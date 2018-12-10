import pytest
from pathlib import Path
from pyWavelet.snowradar import SnowRadar, OIB, AWI

# test the main SnowRadar class by itself plus the OIB and AWI subclasses
TEST_DATA_ROOT = Path(__file__).parent.parent / 'pyWavelet' / 'data'
BASE_CLASS_TEST_FILE = str(TEST_DATA_ROOT / 'sr' / 'Data_20160419_04_010.mat')
OIB_TEST_FILE = str(TEST_DATA_ROOT / 'sr' / 'Data_20160419_04_010.mat')
AWI_TEST_FILE = str(TEST_DATA_ROOT / 'awi' / 'Data_20170410_01_006.mat')

# fixtures are basically objects only instantiated when tests request them specifically
@pytest.fixture
def sr_full():
    return SnowRadar(BASE_CLASS_TEST_FILE, l_case='full')

@pytest.fixture
def sr_meta():
    return SnowRadar(BASE_CLASS_TEST_FILE, l_case='meta')

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

def test_awi_str_repr(awi_full, awi_meta):
    expected = f'AWI Datafile: {Path(AWI_TEST_FILE).name}'
    assert str(awi_full) == expected
    assert str(awi_meta) == expected

#def test_oib_str_repr(oib_full, oib_meta):
def test_oib_str_repr(oib_meta):
    expected = f'OIB Datafile: {Path(OIB_TEST_FILE).name}'
    #assert str(oib_full) == expected
    assert str(oib_meta) == expected

def test_awi_as_dict(awi_full, awi_meta):
    '''expected data based on Data_20170410_01_006.mat'''
    expected_full = {} # placeholder for expected metadata
    expected_meta = {} # placeholder for expected metadata
    '''
    assert awi_full.as_dict() == expected_full
    assert awi_meta.as_dict() == expected_meta
    '''

#def test_oib_as_dict(oib_full, oib_meta):
def test_oib_as_dict(oib_meta):
    '''expected data based on Data_20160419_04_010.mat'''
    expected_full = {} # placeholder for expected metadata
    expected_meta = {} # placeholder for expected metadata
    '''
    assert oib_full.as_dict() == expected_full
    assert oib_meta.as_dict() == expected_meta
    '''
