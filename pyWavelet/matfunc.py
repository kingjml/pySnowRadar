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

def h5todict(h5f, path="/", exclude_names=None):
    '''
    A recursive function to create dictionaries from HDF5/MATLAB7
    '''
    ddict = {}
    for key in h5f[path]:
        if exclude_names is not None:
            if key in exclude_names:
                continue
        if isinstance(h5f[path + "/" + key], h5py._hl.group.Group):
            ddict[key] = h5todict(h5f, path + "/" + key)
        else:
            ddict[key] = np.squeeze(h5f[path + "/" + key][...])
    return ddict