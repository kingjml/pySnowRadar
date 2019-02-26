import pytest
from pathlib import Path
from pySnowRadar.snowradar import SnowRadar

# test the main SnowRadar class by itself plus the OIB and AWI subclasses
TEST_DATA_ROOT = Path(__file__).parent.parent / 'pySnowRadar' / 'data'
OIB_TEST_FILE = TEST_DATA_ROOT / 'sr' / 'Data_20160419_04_010.mat'
AWI_TEST_FILE = TEST_DATA_ROOT / 'awi' / 'Data_20170410_01_006.mat'
NSIDC_TEST_FILE = TEST_DATA_ROOT / 'nsidc' / 'IRSNO1B_20171125_01_118.nc'
TEST_NOT_A_FILE = Path('./fake_sr.mat')

# fixtures are basically objects only instantiated when tests request them specifically
@pytest.fixture
def oib_full():
    '''Full load class for OIB matfile Data_20160419_04_010.mat'''
    return SnowRadar(str(OIB_TEST_FILE), l_case='full')

@pytest.fixture
def oib_meta():
    '''Meta load class for OIB matfile Data_20160419_04_010.mat'''
    return SnowRadar(str(OIB_TEST_FILE), l_case='meta')

@pytest.fixture
def awi_full():
    '''Full load class for AWI matfile Data_20170410_01_006.mat'''
    return SnowRadar(str(AWI_TEST_FILE), l_case='full')

@pytest.fixture
def awi_meta():
    '''Meta load class for AWI matfile Data_20170410_01_006.mat'''
    return SnowRadar(str(AWI_TEST_FILE), l_case='meta')

@pytest.fixture
def nsidc_full():
    '''Full load class for NSIDC matfile IRSNO1B_20171125_01_118.nc'''
    return SnowRadar(str(NSIDC_TEST_FILE), l_case='full')

@pytest.fixture
def nsidc_meta():
    '''Meta load class for NSIDC matfile IRSNO1B_20171125_01_118.nc'''
    return SnowRadar(str(NSIDC_TEST_FILE), l_case='meta')

def test_file_missing():
    with pytest.raises(FileNotFoundError):
        fake_file = SnowRadar('this-file-definitely-does-not-exist-on-anyones-computer-probably.jpeg.gif.tif.ogg.mp4', l_case='full')

def test_bad_file():
    with pytest.raises(IOError):
        bad_file = SnowRadar(TEST_NOT_A_FILE, 'meta')

def test_invalid_load_case():
    with pytest.raises(ValueError):
        wrong_l_case = SnowRadar(str(OIB_TEST_FILE), l_case='not_a_real_l_CASE')

def test_awi_str_repr(awi_full, awi_meta):
    expected = f'AWI_MAT Datafile: {Path(AWI_TEST_FILE).name}'
    assert str(awi_full) == expected
    assert str(awi_meta) == expected

def test_oib_str_repr(oib_full, oib_meta):
    expected = f'OIB_MAT Datafile: {Path(OIB_TEST_FILE).name}'
    assert str(oib_full) == expected
    assert str(oib_meta) == expected
    
def test_nsidc_str_repr(nsidc_full, nsidc_meta):
    expected = f'NSIDC_NC Datafile: {Path(NSIDC_TEST_FILE).name}'
    assert str(nsidc_full) == expected
    assert str(nsidc_meta) == expected

def test_awi_as_dict(awi_full, awi_meta):
    '''expected data based on Data_20170410_01_006.mat'''
    expected_full = {
        'fname': AWI_TEST_FILE.name,
        'fpath': str(AWI_TEST_FILE.absolute()),
        'l_case': 'full',
        'tstart': 1491867198.6980438, # utc time
        'tend': 1491867294.850502, # utc time
        'poly': 'POLYGON ((-156.4123497738317 71.33410559555396, -156.4123497738317 71.3672588206452, -156.503229401445 71.3672588206452, -156.503229401445 71.33410559555396, -156.4123497738317 71.33410559555396))'
    }
    expected_meta = {
        'fname': AWI_TEST_FILE.name,
        'fpath': str(AWI_TEST_FILE.absolute()),
        'l_case': 'meta',
        'tstart': 1491867198.6980438, # utc time
        'tend': 1491867294.850502, # utc time
        'poly': 'POLYGON ((-156.4123497738317 71.33410559555396, -156.4123497738317 71.3672588206452, -156.503229401445 71.3672588206452, -156.503229401445 71.33410559555396, -156.4123497738317 71.33410559555396))'
    }
    assert awi_full.as_dict() == expected_full
    assert awi_meta.as_dict() == expected_meta

def test_oib_as_dict(oib_full, oib_meta):
    '''expected data based on Data_20160419_04_010.mat'''
    expected_full = {
        'fname': OIB_TEST_FILE.name,
        'fpath': str(OIB_TEST_FILE.absolute()),
        'l_case': 'full',
        'tstart': 1461071183.9615262, # utc time
        'tend': 1461071225.002684, # utc time
        'poly': 'POLYGON ((-86.76076421204809 80.22643933420355, -86.76076421204809 80.2711290662714, -86.76363201954442 80.2711290662714, -86.76363201954442 80.22643933420355, -86.76076421204809 80.22643933420355))'
    } 
    expected_meta = {
        'fname': OIB_TEST_FILE.name,
        'fpath': str(OIB_TEST_FILE.absolute()),
        'l_case': 'meta',
        'tstart': 1461071183.9615262, # utc time
        'tend': 1461071225.002684, # utc time
        'poly': 'POLYGON ((-86.76076421204809 80.22643933420355, -86.76076421204809 80.2711290662714, -86.76363201954442 80.2711290662714, -86.76363201954442 80.22643933420355, -86.76076421204809 80.22643933420355))'
    }  
    assert oib_full.as_dict() == expected_full
    assert oib_meta.as_dict() == expected_meta

def test_nsidc_as_dict(nsidc_full, nsidc_meta):
    '''expected data based on IRSNO1B_20171125_01_118.nc'''
    expected_full = {
        'fname': NSIDC_TEST_FILE.name,
        'fpath': str(NSIDC_TEST_FILE.absolute()),
        'l_case': 'full',
        'tstart': 1511659482.0,
        'tend': 1511659516.0,
        'poly': 'POLYGON ((-90.94413070322921 -73.09002326961142, -90.94413070322921 -73.08411595680271, -91.07643835902638 -73.08411595680271, -91.07643835902638 -73.09002326961142, -90.94413070322921 -73.09002326961142))'
    }
    expected_meta = {
        'fname': NSIDC_TEST_FILE.name,
        'fpath': str(NSIDC_TEST_FILE.absolute()),
        'l_case': 'meta',
        'tstart': 1511659482.0,
        'tend': 1511659516.0,
        'poly': 'POLYGON ((-90.94413070322921 -73.09002326961142, -90.94413070322921 -73.08411595680271, -91.07643835902638 -73.08411595680271, -91.07643835902638 -73.09002326961142, -90.94413070322921 -73.09002326961142))'
    }
    assert nsidc_full.as_dict() == expected_full
    assert nsidc_meta.as_dict() == expected_meta

def test_calcpulsewidth_result(awi_full, awi_meta, oib_full, oib_meta, nsidc_full, nsidc_meta):
    awi_full.calcpulsewidth()
    awi_meta.calcpulsewidth()
    # expected AWI values based on Data_20170410_01_006.mat
    expected_awi_null_to_null_pulse_width = 0.075697595645
    expected_awi_equivalent_pulse_width = 0.028389437310606058
    assert awi_full.n2n == expected_awi_null_to_null_pulse_width
    assert awi_full.epw == expected_awi_equivalent_pulse_width
    assert awi_meta.n2n == expected_awi_null_to_null_pulse_width
    assert awi_meta.epw == expected_awi_equivalent_pulse_width
    oib_full.calcpulsewidth()
    oib_meta.calcpulsewidth()
    # expected OIB values based on Data_20160419_04_010.mat
    expected_oib_null_to_null_pulse_width = 0.20186025505333335
    expected_oib_equivalent_pulse_width = 0.07570516616161617
    assert oib_full.n2n == expected_oib_null_to_null_pulse_width
    assert oib_full.epw == expected_oib_equivalent_pulse_width
    assert oib_meta.n2n == expected_oib_null_to_null_pulse_width
    assert oib_meta.epw == expected_oib_equivalent_pulse_width
    # expected NSIDC values based on IRSNO1B_20171125_01_118.nc
    nsidc_full.calcpulsewidth()
    nsidc_meta.calcpulsewidth()
    expected_nsidc_null_to_null_pulse_width = 0.075697595645
    expected_nsidc_equivalent_pulse_width = 0.028389437310606058
    assert nsidc_full.n2n == expected_nsidc_null_to_null_pulse_width
    assert nsidc_full.epw == expected_nsidc_equivalent_pulse_width
    assert nsidc_full.n2n == expected_nsidc_null_to_null_pulse_width
    assert nsidc_meta.epw == expected_nsidc_equivalent_pulse_width
