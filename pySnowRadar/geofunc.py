'''
Functions to handle geospatial manipulation for non-product-specific cases
'''
from pyproj import Proj, transform
from functools import partial

SRC_PROJ = Proj(init='epsg:4326')

def project_latlon(lon, lat, epsg=3413):
    '''
    Project input lon/lat arrays to x/y using the given EPSG code

    Arguments:
        epsg: integer EPSG code for the output coordinate reference system
        (default 3413 for NSIDC North Pole Stereo)

    Outputs:
        2 numpy arrays of projected coordinates: X and Y, respectively
    '''
    try:
        dst_proj = Proj(init=f'epsg:{epsg}')
    except RuntimeError:
        raise ValueError('Supplied EPSG code: %d unrecognized' % epsg)
    trans = partial(transform, SRC_PROJ, dst_proj)
    return trans(lon, lat)
