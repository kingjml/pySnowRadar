import pytest
from pathlib import Path
from pyWavelet.snowradar import SnowRadar, OIB, AWI

# test the main SnowRadar class by itself plus the OIB and AWI subclasses
TEST_DATA_ROOT = Path(__file__).parent.parent / 'pyWavelet' / 'data'
BASE_CLASS_TEST_FILE = TEST_DATA_ROOT / 'sr' / 'Data_20160419_04_010.mat'
OIB_TEST_FILE = TEST_DATA_ROOT / 'sr' / 'Data_20160419_04_010.mat'
AWI_TEST_FILE = TEST_DATA_ROOT / 'awi' / 'Data_20170410_01_006.mat'

# fixtures are basically objects only instantiated when tests request them specifically
@pytest.fixture
def sr_full():
    return SnowRadar(str(BASE_CLASS_TEST_FILE), l_case='full')

@pytest.fixture
def sr_meta():
    return SnowRadar(str(BASE_CLASS_TEST_FILE), l_case='meta')

@pytest.fixture
def oib_full():
    return OIB(str(OIB_TEST_FILE), l_case='full')

@pytest.fixture
def oib_meta():
    return OIB(str(OIB_TEST_FILE), l_case='meta')

@pytest.fixture
def awi_full():
    return AWI(str(AWI_TEST_FILE), l_case='full')

@pytest.fixture
def awi_meta():
    return AWI(str(AWI_TEST_FILE), l_case='meta')

def test_awi_str_repr(awi_full, awi_meta):
    expected = f'AWI Datafile: {Path(AWI_TEST_FILE).name}'
    assert str(awi_full) == expected
    assert str(awi_meta) == expected

def test_oib_str_repr(oib_full, oib_meta):
    expected = f'OIB Datafile: {Path(OIB_TEST_FILE).name}'
    assert str(oib_full) == expected
    assert str(oib_meta) == expected

def test_base_as_dict(sr_full, sr_meta):
    expected_full = {
        'fname': BASE_CLASS_TEST_FILE.name,
        'fpath': str(BASE_CLASS_TEST_FILE.absolute()),
        'l_case': 'full'
    }
    expected_meta = {
        'fname': BASE_CLASS_TEST_FILE.name,
        'fpath': str(BASE_CLASS_TEST_FILE.absolute()),
        'l_case': 'meta'
    }
    assert sr_full.as_dict() == expected_full
    assert sr_meta.as_dict() == expected_meta

def test_awi_as_dict(awi_full, awi_meta):
    '''expected data based on Data_20170410_01_006.mat'''
    expected_full = {
        'fname': AWI_TEST_FILE.name,
        'fpath': str(AWI_TEST_FILE.absolute()),
        'l_case': 'full',
        'tstart': 1491867198.6980438,
        'tend': 1491867294.850502,
        'poly': 'POLYGON ((-156.4123497738317 71.33410559555396, -156.4123497738317 71.3672588206452, -156.503229401445 71.3672588206452, -156.503229401445 71.33410559555396, -156.4123497738317 71.33410559555396))'
    }
    expected_meta = {
        'fname': AWI_TEST_FILE.name,
        'fpath': str(AWI_TEST_FILE.absolute()),
        'l_case': 'meta',
        'tstart': 1491867198.6980438,
        'tend': 1491867294.850502,
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
        'tstart': 1461071183.9615262,
        'tend': 1461071225.002684,
        'poly': 'POLYGON ((-86.76076421204809 80.22643933420355, -86.76076421204809 80.2711290662714, -86.76363201954442 80.2711290662714, -86.76363201954442 80.22643933420355, -86.76076421204809 80.22643933420355))'
    } 
    expected_meta = {
        'fname': OIB_TEST_FILE.name,
        'fpath': str(OIB_TEST_FILE.absolute()),
        'l_case': 'meta',
        'tstart': 1461071183.9615262,
        'tend': 1461071225.002684,
        'poly': 'POLYGON ((-86.76076421204809 80.22643933420355, -86.76076421204809 80.2711290662714, -86.76363201954442 80.2711290662714, -86.76363201954442 80.22643933420355, -86.76076421204809 80.22643933420355))'
    }  
    assert oib_full.as_dict() == expected_full
    assert oib_meta.as_dict() == expected_meta

def test_calcpulsewidth_result(awi_full, awi_meta, oib_full, oib_meta):
    '''
    don't even THINK about testing the SnowRadar parent class, because 
    it doesn't have a self.bandwidth attribute
    '''
    awi_full.calcpulsewidth()
    awi_meta.calcpulsewidth()
    oib_full.calcpulsewidth()
    oib_meta.calcpulsewidth()