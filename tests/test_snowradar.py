import pytest
from pathlib import Path
from pyWavelet.snowradar import OIB, AWI

# test the main class by itself plus the OIB and AWI subclasses

TEST_DATA_ROOT = Path(__file__).parent.parent / 'pyWavelet' / 'data'
OIB_TEST_FILE = str(TEST_DATA_ROOT / 'sr' / 'Data_20160419_04_010.mat')
AWI_TEST_FILE = str(TEST_DATA_ROOT / 'awi' / 'Data_20170410_01_006.mat')

# these are the testing instances that will be used for further test functions
@pytest.fixture
def oib_full():
    return OIB(OIB_TEST_FILE, l_case='full')

@pytest.fixture
def oib_meta():
    return OIB(OIB_TEST_FILE, l_case='meta')

@pytest.fixture
def awi_full():
    return AWI(AWI_TEST_FILE, l_case='full')

@pytest.fixture
def awi_meta():
    return AWI(AWI_TEST_FILE, l_case='meta')

def test_awi_full_thing(oib_full):
    assert type(str(awi_full)) == str

def test_snowradar_as_dict(oib_full, oib_meta):
    print(oib_meta.time_gps)
    assert oib_full.as_dict() == None
    '''
    assert oib_meta.as_dict() == {
        'fname': oib_meta.file_name,
        'fpath': oib_meta.file_path,
        'tstart': oib_meta.time_utc[0],
        'tend': oib_meta.time_utc[-1],
        'poly': oib_meta.poly
    }
    '''