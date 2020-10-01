'''
Functions to handle geospatial manipulation for non-product-specific cases
'''
from pyproj import Proj, transform
from functools import partial
import numpy as np

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

def haversine(lat, lon, to_radians=True, earth_radius=6371):
    if to_radians:
        lat1, lon1 = np.radians([lat, lon])
        lat2 = np.roll(lat1,1)
        lon2 = np.roll(lon1,1)

    a = np.sin((lat1-lat2)/2.0)**2 + \
        np.cos(lat2) * np.cos(lat1) * np.sin((lon1-lon2)/2.0)**2
    
    a[0] = np.nan
    
    return earth_radius * 2 * np.arcsin(np.sqrt(a))