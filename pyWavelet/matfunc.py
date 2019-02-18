'''
Functions to handle MATLAB data formats
MAT file reader functions from here: https://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries
'''
import scipy.io as spio
import h5py
import numpy as np

def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
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
    A recursive function which constructs nested dictionaries from matobjects 
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
    '''Loads SnowRadar MAT files depending on source'''
    # The versions of the OIB MAT file formats seem to be inconsistent
    # This tries the two different approaches to read to find one that works
    try:
        radar_dat = loadmat(sr_obj.file_path)
    except NotImplementedError:
        with h5py.File(sr_obj.file_path, 'r') as in_h5:
            radar_dat = h5py_to_dict(in_h5, exclude_names='#refs#')
    except:
        raise IOError('Could not read SnowRadar file: %s' % sr_obj.file_name)
    return radar_dat


def h5py_to_dict(hdf5_obj, exclude_names=None):
    '''
    For the HDF5-compliant SnowRadar datasets (OIB 2017, AWI), this reader method
    works for the majority of the important data.

    There are still some HDF5-Object-References that do not get parsed for whatever
    reason
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
                if v.value.dtype == 'uint16':
                    # converts uint16 to Byte and decodes to string
                    # https://stackoverflow.com/a/45593385
                    data[k] = v.value.tobytes()[::2].decode()
                elif v.value.dtype == 'uint8':
                    # converts 1x1 uint8 array to single bool value
                    data[k] = v.value.astype(bool).item()
                # empty attribute stored as 2-element array of type uint64
                elif v.value.dtype == 'uint64' and v.value.shape == (2,):
                    data[k] == None
            else:
                # need to convert HDF object references into actual data
                dtype = v.value.dtype
                if dtype == 'object':
                    try:
                        data[k] = np.array([
                            np.squeeze(v.parent[reference].value)
                            for reference in np.squeeze(v.value)
                        ])
                    except:
                        pass
                else:
                    try: # try to convert 1x1 arrays into scalars
                        data[k] = v.value.item()
                    except ValueError: # keep NxN arrays the way they are
                        data[k] = np.squeeze(v.value)
        elif isinstance(v, h5py.Group):
            # recursively step through h5py.Groups to get at the datasets
            data[k] = h5py_to_dict(v)
    return data