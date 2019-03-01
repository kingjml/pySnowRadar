from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from pySnowRadar.picklayers import picklayers
from pySnowRadar.snowradar import SnowRadar


def geo_filter(input_sr_data):
    '''
    Given a list of SnowRadar datafiles (.mat, .h5, .nc), filter out
    any files whose bounding geometry intersect with the Canadian land mass

    Landmask is based on Natural Earth 110m (low-res)

    Arguments:
        input_sr_data: list of supported SnowRadar data files

    Output:
        subset of input_sr_data where no land intersections occur
    '''
    # Drop all data that intersects with Canadian land features
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    canada = world.loc[world.name == 'Canada']
    # Load the datafiles in 'meta' mode to just scrape the bounding geometry
    sr_meta = [SnowRadar(sr, 'meta') for sr in input_sr_data]
    sr_gdf = gpd.GeoDataFrame(
        data={'file': [sr.file_path for sr in sr_meta]}, 
        geometry=[sr.poly for sr in sr_meta], 
        crs={'init':'epsg:4326'}
    )
    sr_gdf = sr_gdf.drop(
        gpd.sjoin(sr_gdf, canada, how='inner', op='intersects').index
    )
    if len(sr_gdf) < 1:
        print('No suitable datafiles left after geospatial filtering')
        return []
    return sr_gdf.file.tolist()


def extract_layers(data_path, snow_density=0.3, dump_results=False):
    '''
    For a given SnowRadar datafile, estimate the air-snow and snow-ice interfaces
    based on the supplied snow density

    Arguments:
        data_path: file path to input SnowRadar data file
        snow_density: a single float value to use when picking interface layers
        dump_results: whether or not to save the dataframe to a local csv file under ./dump/

    Output:
        A pandas dataframe with the following columns:
            'src': the name of the source SnowRadar data file
            'lat': latitude of trace(?)
            'lon': longitude of trace(?)
            'n_snow': the refractive index used during layer picking
            'b_as': the picked air-snow interface layer
            'b_si': the picked snow-ice interface layer
            'snow_depth': estimated snow depth based on picked layers
    '''
    # skip reprocessing if local csv dump exists
    outpath = Path('./dump')
    outname = Path(data_path).stem + '.csv'
    outfile = outpath / outname
    if outfile.exists():
        print('File exists for %s. Skipping processing....' % Path(data_path).name)
        result = pd.read_csv(str(outfile.stem), index_col=0)
        return result

    if not(0.1 <= snow_density <= 0.4):
        raise ValueError('Invalid snow density passed: %.3f ' % snow_density + \
                        '(Must be between 0.1 and 0.4)')
    # Convert density to refractive index
    r_idx = np.sqrt((1 + 0.51 * snow_density) ** 3)
    radar_dat = SnowRadar(data_path, 'full')
    radar_dat.surf_bin, radar_dat.surface = radar_dat.get_surface()
    radar_dat.calcpulsewidth()
    lower, upper = radar_dat.get_bounds(m_above=5)
    radar_subset = radar_dat.data_radar[upper:lower, :]
    try:
        airsnow, snowice = np.apply_along_axis(
            picklayers,
            axis=0, 
            arr=radar_subset,
            null_2_space=radar_dat.n2n,
            delta_fast_time_range=radar_dat.dfr,
            n_snow=r_idx
        )
    except:
        # Set interfaces to NaN if anything goes wrong
        print('Picklayers blew up for file: %s' % radar_dat.file_name)
        airsnow = np.array([np.nan] * radar_dat.lat.shape[0])
        snowice = np.array([np.nan] * radar_dat.lat.shape[0])

    # Some extra columns
    snow_depth = (snowice - airsnow) * radar_dat.dfr / r_idx
    data_src = np.array([radar_dat.file_name] * radar_dat.lat.shape[0])
    refractive_index = np.array([r_idx] * radar_dat.lat.shape[0])

    result = pd.DataFrame({
        'src': data_src,
        'lat': radar_dat.lat,
        'lon': radar_dat.lon,
        'n_snow': refractive_index,
        'b_as': airsnow,
        'b_si': snowice,
        'snow_depth': snow_depth
    })

    if dump_results:
        outpath.mkdir(parents=True, exist_ok=True)
        result.to_csv(str(outfile), na_rep='nan')
    
    return result




def batch_process(input_sr_data, snow_density=0.3, workers=4, dump_results=True):
    '''
    For a given list of SnowRadar data file paths:
        1) Pick air-snow and snow-ice interface layers for the 
           data files using the supplied snow density

        2) Produce a dataframe with the all the picked interfaces for each 
           of the data files

    Arguments:
        input_sr_data: list of supported SnowRadar data files
        snow_density: 
            a) a single float value to apply to all data files
            b) a list of file-specific density values matching the length of input_sr_data
        workers: number of worker processes to use
        dump_results: dumps each dataframe to a local csv

    Output:
        A concatendated pandas dataframe with the following columns:
            'src': the name of the source SnowRadar data file
            'lat': latitude of trace(?)
            'lon': longitude of trace(?)
            'n_snow': the refractive index used during layer picking
            'b_as': the picked air-snow interface layer
            'b_si': the picked snow-ice interface layer
            'snow_depth': estimated snow depth based on picked layers

    '''
    # first check to ensure that you're not going to blow up your hardware
    current_cores = cpu_count()
    if workers > current_cores:
        raise SystemError('workers argument (passed: %d) cannot ' % workers + 
                          'exceed current CPU count (%d)' % current_cores)
    length = len(input_sr_data)
    dump_triggers = [dump_results] * length                          
    # Single value passed for snow_density
    if isinstance(snow_density, float):
        if not(0.1 <= snow_density <= 0.4):
            raise ValueError('Invalid snow density passed: %.3f ' % snow_density + \
                            '(Must be between 0.1 and 0.4)')
        else:
            process_args = zip(
                input_sr_data, 
                [snow_density] * length, 
                dump_triggers
            )
    # Multiple values passed for snow_density
    elif isinstance(snow_density, list):
        if len(snow_density) != length:
            raise ValueError('Passed list of snow densities must match length of ' + \
                             'passed snowradar files. ' + \
                             'Files -> %d Densities -> %d' % (len(snow_density), length))
        elif any([not(0.1 <= sd <= 0.4) for sd in snow_density]):
            raise ValueError('Invalid list of snow densities passed. ' + \
                             'All snow density values must be between 0.1 and 0.4')
        else:
            process_args = zip(
                input_sr_data,
                snow_density,
                dump_triggers
            )
    else:
        print('Invalid argument for snow density. Proceeding with default value (0.3)')
        process_args = zip(
            input_sr_data, 
            [snow_density] * length,
            dump_triggers
        )

    # Define the threadpool and submit/execute the list of tasks
    with ProcessPoolExecutor(workers) as pool:
        futures = [pool.submit(extract_layers, *foo) for foo in process_args]
        results = [f.result() for f in futures]
    # return a concatenated dataframe containing results for all input datasets
    df = pd.concat(results)
    return df
