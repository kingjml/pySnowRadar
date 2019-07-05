'''
Functions to handle MATLAB and NetCDF data formats
MAT file reader functions from here: https://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries
'''
import scipy.io as spio
import h5py
from dateutil.parser import parse
import numpy as np

def loadmat(filename):
    '''
    this custom function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to convert all entries
    which are still mat-objects
    '''
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(in_dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in in_dict:
        if isinstance(in_dict[key], spio.matlab.mio5_params.mat_struct):
            in_dict[key] = _todict(in_dict[key])
    return in_dict        

def _todict(matobj):
    '''
    A recursive function which constructs nested dictionaries from mat-objects 
    '''
    out_dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            out_dict[strg] = _todict(elem)
        else:
            out_dict[strg] = elem
    return out_dict


def unified_loader(sr_obj):
    '''
    The entry-point called by the SnowRadar class, that determines 
    how each different type of SnowRadar data file is parsed
    
    Currently operates in the following decision tree:

    Does the data filename end in '.nc'?
        Yes: Treat as a NSIDC OIB L1b NetCDF    - load using h5py
        No: Can the scipy loadmat() function process the data file?
            Yes: Treat as MatLab v5             - load using scipy 
            No: Treat as MatLab v7 / AWI        - load using h5py

    Arguments:
        sr_obj: an instance of SnowRadar class

    Output:
        radar_dat: a dictionary containing the data extracted from 
                   the input SnowRadar instance's source dataset (sr_obj.file_path)
        
    '''
    # Check for NSIDC NetCDF format first
    if sr_obj.file_path.endswith('.nc'):
        try:
            with h5py.File(sr_obj.file_path, 'r') as in_nc:
                radar_dat = nc_to_dict(in_nc)
        except:
            raise IOError('Could not read SnowRadar file: %s' % sr_obj.file_name)
    # else assume that it can only be a Matlab data file
    else:
        # The versions of the OIB MAT file formats seem to be inconsistent
        # This tries different approaches to read to find one that works
        try:
            radar_dat = loadmat(sr_obj.file_path)
        # NotImplementedError raised by spio.loadmat for matlab v7 files
        except NotImplementedError:
            with h5py.File(sr_obj.file_path, 'r') as in_mat:
                radar_dat = h5py_to_dict(in_mat, exclude_names='#refs#')
        except:
            raise IOError('Could not read SnowRadar file: %s' % sr_obj.file_name)
    return radar_dat


def h5py_to_dict(hdf5_obj, exclude_names=None):
    '''
    For the MatLab-but-HDF5-compliant SnowRadar datasets (OIB 2017, AWI), 
    this reader method works for the majority of the important data.

    There are still some HDF5-Object-References that do not get parsed 
    for unknown reasons (possibly due to not catching certain variable-types)

    Arguments:
        hdf5_obj: an open h5py.File (preferably in read-mode)
        exclude_names: a string or list of strings of specific top-level
                       dataset categories to ignore

    Output:
        data: a dictionary containing the data extracted from hdf5_obj
    '''
    data = {}
    # Note that exclude_names can be tricky since it
    # is applied recursively when it hits a h5py.Group
    if not type(exclude_names) == list:
        exclude_names = [exclude_names]
    for k, v in hdf5_obj.items():
        if k in exclude_names:
            continue
        if isinstance(v, h5py.Dataset):
            # some important strings are encoded as uint16 for whatever reason
            # thankfully there is a unique HDF attribute that flags encoded values
            #
            # Note: this doesn't decode all the strings, but it does work for the 
            # data we need
            must_decode = 'MATLAB_int_decode' in list(v.attrs)
            if must_decode:
                # treat uint8 differently from uint16
                if v[()].dtype == 'uint16':
                    # converts uint16 to Byte and decodes to string
                    # https://stackoverflow.com/a/45593385
                    data[k] = v[()].tobytes()[::2].decode()
                elif v[()].dtype == 'uint8':
                    # converts 1x1 uint8 array to single bool value
                    data[k] = v[()].astype(bool).item()
                # empty attribute stored as 2-element array of type uint64
                elif v[()].dtype == 'uint64' and v[()].shape == (2,):
                    data[k] == None
            else:
                # need to convert HDF object references into actual data
                dtype = v[()].dtype
                if dtype == 'object':
                    try:
                        data[k] = np.array([
                            np.squeeze(v.parent[reference][()])
                            for reference in np.squeeze(v[()])
                        ])
                    except:
                        pass
                else:
                    try: # try to convert 1x1 arrays into scalars
                        data[k] = v[()].item()
                    except ValueError: # keep NxN arrays the way they are
                        data[k] = np.squeeze(v[()])
        elif isinstance(v, h5py.Group):
            # recursively step through h5py.Groups to get at the datasets
            data[k] = h5py_to_dict(v)
    return data

# because the NC structure doesn't match the Matlab structure
# and we've already written all the code to match the Matlab structure...
NC_VAR_NAME_SWAP_LUT = {
    'amplitude': 'Data',
    'lat': 'Latitude',
    'lon': 'Longitude',
    'roll': 'Roll',
    'pitch': 'Pitch',
    'fasttime': 'Time',
    'altitude': 'Elevation', # the NSIDC user guide lists elevation as 'alt'
    'alt': 'Elevation',      # but the downloaded NC file uses 'altitude'
    'time': 'UTC_Time'
}

def nc_to_dict(hdf5_obj):
    '''
    Using h5py, extract a subset of the input HDF dataset's 
    attributes and data arrays

    NB: this is currently hardcoded to the current (Feb 26 2019)
    NSIDC OIB L1b NetCDF datasets

    Arguments:
        hdf5_obj: an open h5py.File (preferably in read-mode)

    Output:
        data: a dictionary containing the data extracted from hdf5_obj
    '''
    data = {}
    # iterate through the netCDF, trying to convert to dict
    # in the same format as the Matlab dicts
    for k, v in hdf5_obj.items():
        k = NC_VAR_NAME_SWAP_LUT.get(k, k)
        k_list = [s.split('(')[0].split('{')[0] for s in k.split('.')]
        if v.dtype == '<S1':
            try:
                val = ''.join(v[()].flatten().astype('U'))
            except TypeError:
                val = None
        elif v.dtype == 'uint16':
            val = v[()].astype(bool).item()
        else:
            val = np.squeeze(v[()])
        nc_set_nested_value(data, k_list, val)
    # also grab the GPS times from the netCDF header
    min_time = parse(
        hdf5_obj.attrs['min_time_bound'].astype('U').rstrip(' GPS')
    ).timestamp()
    max_time = parse(
        hdf5_obj.attrs['max_time_bound'].astype('U').rstrip(' GPS')
    ).timestamp()
    data['GPS_time'] = np.array([min_time, max_time], dtype=float)
    # Altering the NSIDC L1b datasets to conform to Matlab data
    data['Data'] = np.power(10, data['Data'][:] / 10)
    data['Time'] *= 1e-6
    data['Roll'] = np.radians(data['Roll'][:])
    data['Pitch'] = np.radians(data['Pitch'][:])
    return data

def nc_set_nested_value(target_dict, keys, value):
    '''
    Given a target dictionary, insert a nested dictionary 
    keyed to the list of keys, with the deepest key holding
    the passed value

    https://stackoverflow.com/a/21298015
    '''
    deepest_key = keys.pop()
    for key in keys:
        target_dict = target_dict.setdefault(key, {})
    target_dict.setdefault(deepest_key, value)
