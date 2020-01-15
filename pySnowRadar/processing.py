import os
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
#from pySnowRadar.picklayers import picklayers #obsolete picklayers function
from pySnowRadar import ATM, SnowRadar, pick_layers


def geo_filter(input_sr_data):
    '''
    Given a list of SnowRadar datafiles (.mat, .h5, .nc), filter out
    any files whose bounding geometry intersect with land

    Landmask is based on NaturalEarth 1:10m Cultural v4.1.0 (Canada, Greenland, and USA)
    http://www.naturalearthdata.com/

    Arguments:
        input_sr_data: list of supported SnowRadar data files

    Output:
        subset of input_sr_data where no land intersections occur
    '''
    # Drop all data that intersects with land features
    land = gpd.read_file('/vsizip/' + str(Path(__file__).parent / 
                         'data' / 'natearth' / 
                         'ne_10m_admin_0_countries_northamerica.zip'))
    # Load the datafiles in 'meta' mode to just scrape the simplified track line
    sr_meta = [SnowRadar(sr, 'meta') for sr in input_sr_data]
    sr_gdf = gpd.GeoDataFrame(
        data={'file': [sr.file_path for sr in sr_meta]}, 
        geometry=[sr.line for sr in sr_meta], 
        crs={'init':'epsg:4326'}
    )
    sr_gdf = sr_gdf.drop(
        gpd.sjoin(sr_gdf, land, how='inner', op='intersects').index
    )
    if len(sr_gdf) < 1:
        print('No suitable datafiles left after geospatial filtering')
        return []
    return sr_gdf.file.tolist()


def extract_layers(data_path, snow_density=0.3, picker='Wavelet-TN', dump_results=False):
    '''
    For a given SnowRadar datafile, estimate the air-snow and snow-ice interfaces
    using the supplied picker and snow density

    NB: This is currently only working for the Wavelet-TN picker method. There will probably
    need to be restructuring for any other picker algorithm (i.e. one that requires more than
    just n2n, dfr, n_snow, and snow density)

    Arguments:
        data_path: file path to input SnowRadar data file
        snow_density: a single float value to use when picking interface layers
        picker: which picker algorithm to apply (Wavelet-TN, NSIDC, or GSFC-NK)
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
    if dump_results:
        outpath = Path('./dump')
        outname = Path(data_path).stem + '.csv'
        outfile = outpath / outname
        if outfile.exists():
            print('File exists for %s. Skipping processing....' % Path(data_path).name)
            result = pd.read_csv(str(outfile), index_col=0)
            return result

    if not(0.1 <= snow_density <= 0.4):
        raise ValueError(
            'Invalid snow density passed: %.3f (Must be between 0.1 and 0.4)' % snow_density
        )
    # Convert density to refractive index
    r_idx = np.sqrt((1 + 0.51 * snow_density) ** 3)
    radar_dat = SnowRadar(data_path, 'full')
    radar_dat.surf_bin, radar_dat.surface = radar_dat.get_surface()
    radar_dat.calcpulsewidth()
    lower, upper = radar_dat.get_bounds(m_above=5)
    radar_subset = radar_dat.data_radar[upper:lower, :]
    try:
        # this may be replaceable with a ThreadPoolExecutor... 
        airsnow, snowice = np.apply_along_axis(
            pick_layers,
            axis=0, 
            arr=radar_subset,
            params={
                'n2n':radar_dat.n2n,
                'dfr': radar_dat.dfr,
                'n_snow': r_idx
            },
            picker_func=picker
        )
    except:
        # Set interfaces to NaN if anything goes wrong
        print('pick_layers blew up for file: %s' % radar_dat.file_name)
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


def batch_process(input_sr_data, snow_density=0.3, picker='Wavelet-TN', workers=4, dump_results=True):
    '''
    For a given list of SnowRadar data file paths:
        1) Pick air-snow and snow-ice interface layers for the 
           data files using the provided picker and snow density

        2) Produce a dataframe with the all the picked interfaces for each 
           of the data files

    Arguments:
        input_sr_data: list of supported SnowRadar data files
        snow_density: 
            a) a single float value to apply to all data files
            b) a list of file-specific density values matching the length of input_sr_data
        picker: which picker algorithm to apply (Wavelet-TN, NSIDC, or GSFC-NK)
        workers: number of worker processes to use
        dump_results: dumps each dataframe to a local csv

    Output:
        A concatenated pandas dataframe with the following columns:
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
    picker_args = [picker] * length                     
    # Single value passed for snow_density
    if isinstance(snow_density, float):
        if not(0.1 <= snow_density <= 0.4):
            raise ValueError('Invalid snow density passed: %.3f ' % snow_density + \
                            '(Must be between 0.1 and 0.4)')
        else:
            process_args = zip(
                input_sr_data, 
                [snow_density] * length,
                picker_args, 
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
                picker_args,
                dump_triggers
            )
    else:
        print('Invalid argument for snow density. Proceeding with default value (0.3)')
        process_args = zip(
            input_sr_data, 
            [0.3] * length,
            picker_args,
            dump_triggers
        )

    # Define the threadpool and submit/execute the list of tasks
    with ProcessPoolExecutor(workers) as pool:
        futures = [pool.submit(extract_layers, *foo) for foo in process_args]
        results = [f.result() for f in futures]
    # return a concatenated dataframe containing results for all input datasets
    df = pd.concat(results)
    return df

def fetch_atm_data(sr, atm_folder):
    '''
    Attempt to find and load any locally-available NASA ATM data granules 
    that share the same day as the passed SnowRadar object

    Inputs:
        atm_folder: the local directory where ATM granules should be found
    
    Outputs:
        A dataframe containing concatenated ATM data (if multiple local ATM files exist)
    '''
    if not os.path.isdir(atm_folder):
        raise FileNotFoundError('Cannot locate ATM folder: %s' % os.path.abspath(atm_folder))

    # check for temporal match (same day as current SnowRadar data)
    d = sr.day.strftime('%Y%m%d')
    relevant_atm_data = [
        ATM(os.path.join(r, f)) 
        for r, ds, fs in os.walk(atm_folder) 
        for f in fs if 
        'ATM' in f and 
        f.endswith('.h5') and 
        f.split('_')[1] == d
    ]
    if len(relevant_atm_data) == 0:
        print('No ATM data found for %s' % str(sr))
        return
    
    # check for spatial match (very rough due to simplicity of atm.bbox)
    relevant_atm_data = [
        atm for atm in relevant_atm_data
        if atm.bbox.intersects(sr.line)
    ]
    if len(relevant_atm_data) == 0:
        print('No ATM data found for %s' % str(sr))
        return

    # assuming we still have some ATM data after spatiotemporal filtering, 
    # we sort by filename and concatenate into one big dataframe
    relevant_atm_data.sort(key=lambda x: x.file_name)
    df = pd.concat([
        pd.DataFrame({
            'atm_src': [atm.file_name]*len(atm.pitch),
            'atm_lat': atm.latitude,
            'atm_lon': atm.longitude,
            'atm_elev': atm.elevation,
            'atm_pitch': atm.pitch,
            'atm_roll': atm.roll,
            'atm_time_gps': atm.time_gps
        })
        for atm in relevant_atm_data
    ]).reset_index(drop=True)
    return df
