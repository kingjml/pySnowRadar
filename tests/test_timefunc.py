import pytest
from pyWavelet.timefunc import utcleap

def test_utcleap_invalid():
    with pytest.raises(Exception):
        result = utcleap('a')

def test_utcleap_valid():
    true_time = 1092121230.0
    assert utcleap([1092121243.0]) == true_time