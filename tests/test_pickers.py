import numpy as np
from pathlib import Path
from pySnowRadar import SnowRadar
from pySnowRadar.algorithms import Wavelet_TN, GSFC_NK, NSIDC

# adapted from Data_20160419_04_010.mat
sr = SnowRadar(Path('__file__').parent.parent / 'pySnowRadar' / 'data' / 'sr' / 'Data_20160419_04_010.mat', l_case='full')
sr.surf_bin, sr.surface = sr.get_surface()
bnds = sr.get_bounds(5)
# only use the first trace
data = sr.data_radar[bnds[1]:bnds[0], 0]
n2n = 0.20186025505333335
dfr = 0.012975781596299215
# density set to 0.3 kg/m^3
n_snow = np.sqrt((1 + 0.51 * 0.300) ** 3) 

def test_wavelet_tn():
    results = Wavelet_TN(data, n2n, dfr, n_snow, 1, 10)
    assert len(results) == 2
    airsnow, snowice = results
    assert type(airsnow) == np.int64
    assert type(snowice) == np.int64
    assert snowice > airsnow or snowice == airsnow

def test_gsfc_nk():
    results = GSFC_NK(data)
    assert len(results) == 2
    airsnow, snowice = results
    assert type(airsnow) == np.int64
    assert type(snowice) == np.int64
    assert snowice > airsnow or snowice == airsnow

def test_nsidc():
    results = NSIDC(data)
    assert len(results) == 2
    airsnow, snowice = results
    assert type(airsnow) == np.int64
    assert type(snowice) == np.int64
    assert snowice > airsnow or snowice == airsnow