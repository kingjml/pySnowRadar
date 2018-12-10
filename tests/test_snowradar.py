import pytest
from pathlib import Path
from pyWavelet.snowradar import OIB, AWI

# test the main class by itself plus the OIB and AWI subclasses

TEST_DATA_ROOT = Path(__file__).parent.parent / 'pyWavelet' / 'data'
OIB_TEST_FILE = str(TEST_DATA_ROOT / 'sr' / 'Data_20160419_04_010.mat')
AWI_TEST_FILE = str(TEST_DATA_ROOT / 'awi' / 'Data_20170410_01_006.mat')

# fixtures are basically instantiated objects for testing purposes
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

def test_awi_str_repr(awi_full):
    expected = f'AWI Datafile: {Path(AWI_TEST_FILE).name}'
    assert str(awi_full) == expected

def test_oib_str_repr(oib_full):
    expected = f'OIB Datafile: {Path(OIB_TEST_FILE).name}'
    assert str(oib_full) == expected

def test_awi_as_dict(awi_full, awi_meta):
    assert awi_full.as_dict() == None
    '''
    assert awi_meta.as_dict() == {
        'fname': awi_meta.file_name,
        'fpath': awi_meta.file_path,
        'tstart': awi_meta.time_utc[0],
        'tend': awi_meta.time_utc[-1],
        'poly': awi_meta.poly
    }
    '''

def test_oib_as_dict(oib_full, oib_meta):
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
