import os
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from pySnowRadar import ATM, SnowRadar, algorithms

def calc_snr(sr, noise_bins = 100, interfaces = None):
    '''
    Calculate the signal-to-noise ratio for an entire frame or at 
    specific interfaces

    Input:
        noise_bins: The number of leading bins to consider as noise
        interfaces: Optional numpy array of interface bins to return SNR
    
    Output:
        SNR for an entire frame or at specific bins (interfaces)
    '''
    noise = np.median(sr.data_radar[0: noise_bins, :], axis = 0)
    if interfaces is not None:
        return 10*np.log10((np.array([sr.data_radar[idx, trace] \
                           for trace,idx in enumerate(interfaces)]) / noise))
    else:
        return 10*np.log10(sr.data_radar/noise)

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

def extract_layers(data_path, picker=algorithms.Wavelet_TN, params=None, dump_results=False):
    '''
    For a given SnowRadar datafile, estimate the air-snow and snow-ice interfaces
    using the supplied picker and snow density

    Arguments:
        data_path: file path to input SnowRadar data file
        picker: Which picker algorithm to apply (default is algorithms.Wavelet_TN)
        params: A dictorary of expected parameters to pass to the picker
        dump_results: whether or not to save the dataframe to a local csv file under ./dump/

    Output:
        A pandas dataframe with the following columns:
            'src': the name of the source SnowRadar data file
            'picker': the name of the picker algo
            'lat': latitude of trace
            'lon': longitude of trace
            'n_snow': the refractive index used during layer picking
            'b_ref': the reference bin considered as 0 in the origional file
            'b_as': the picked air-snow interface layer
            'b_si': the picked snow-ice interface layer
            'snow_depth': estimated snow depth based on picked layers
            'params': a dict of all input and generated params
    '''
    
    if dump_results: 
        outpath = Path('./dump')
        outname = Path(data_path).stem + '.csv'
        outfile = outpath / outname
        if outfile.exists(): # skip reprocessing if local csv dump exists
            print('File exists for %s. Skipping processing....' % Path(data_path).name)
            result = pd.read_csv(str(outfile), index_col=0)
            return result
    
    # Check that the picker passed exists
    if not(picker in algorithms.available_pickers()):
        raise ValueError(
            'Invalid picker name:' % picker.__name__
        )

    # TODO: Modify to allow refractive index (n_snow) as an alternative
    if 'snow_density' not in params:
        raise ValueError(
            'Snow density or refractive index input required for all pickers'
        )
    elif (not(0.1 <= params['snow_density'] <= 0.4)):
          raise ValueError(
            'Invalid snow density passed: %.3f (Must be between 0.1 and 0.4)' % snow_density)
    
    # Load radar data 
    radar_dat = SnowRadar(data_path, 'full')
    radar_dat.surf_bin, radar_dat.surface = radar_dat.get_surface()
    radar_dat.calcpulsewidth()
    
    # Subset radar traces to reduce computational load
    # TODO: Should we allow the subset bounds to be user defined?
    lower, upper = radar_dat.get_bounds(m_above=5)
    radar_sub = radar_dat.data_radar[upper:lower, :]
   
    # Calc or init other necessary params
    params['n_snow'] = np.sqrt((1 + 0.51 * params['snow_density']) ** 3)
    params['null_2_space']  = radar_dat.n2n
    params['delta_fast_time_range'] = radar_dat.dfr
    
    # Apply the picker to the file, trace by trace
    try:
        airsnow, snowice = np.apply_along_axis(
            picker, 
            axis=0, 
            arr=radar_sub, 
            **params
        )
    except:
         # Set interfaces to NaN if anything goes wrong
         # TODO: We should catch and print the exception
         print('pick_layers blew up for file: %s' % radar_dat.file_name)
         airsnow = np.array([np.nan] * radar_dat.lat.shape[0])
         snowice = np.array([np.nan] * radar_dat.lat.shape[0])
        
    # Calc snow depth from pics and add the 
    snow_depth = (snowice - airsnow) * radar_dat.dfr / params['n_snow']
    data_src = np.array([radar_dat.file_name] * radar_dat.lat.shape[0])

    result = pd.DataFrame({
        'src': data_src,
        'picker': picker.__name__,
        'lat': radar_dat.lat,
        'lon': radar_dat.lon,
        'n_snow': params['n_snow'],
        #'b_ref': upper,
        'b_as': upper + airsnow,
        'b_si': upper + snowice,
        'snow_depth': snow_depth,
    })
    
    if dump_results:
        outpath.mkdir(parents=True, exist_ok=True)
        result.to_csv(str(outfile), na_rep='nan')
    
    return result

def batch_process(input_sr_data, picker, params, workers=4, dump_results=False):
    '''
    For a given list of SnowRadar data file paths:
        1) Pick air-snow and snow-ice interface layers for the 
           data files using the provided picker and snow density

        2) Produce a dataframe with the all the picked interfaces for each 
           of the data files

    Arguments:
        picker: which picker algorithm to apply (Wavelet-TN, NSIDC, or GSFC-NK)
        params: dictonary of picker paramters or list of dicts
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
    current_cores = cpu_count()
    if workers > current_cores:
        raise SystemError('workers argument (passed: %d) cannot ' % workers + 
                          'exceed current CPU count (%d)' % current_cores)

    # Generate the picker input and output paramters
    length = len(input_sr_data)
    dump_triggers = [dump_results] * length     
    picker_args = [picker] * length                     
    
    if isinstance(params, dict):
    # If the input parameters the same for every file
        process_args = zip(
            input_sr_data, 
            [picker] * length,
            [params] * length,
            dump_triggers
        )
    # If the input parameters vary
    elif isinstance(params, list):
        process_args = zip(
            input_sr_data, 
            [picker] * length,
            params,
            dump_triggers)
    
    with ProcessPoolExecutor(workers) as pool:
        futures = [pool.submit(extract_layers, *foo) for foo in process_args]
        results = [f.result() for f in futures]
    
    # return a concatenated dataframe containing results for all input datasets
    return pd.concat(results)

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
