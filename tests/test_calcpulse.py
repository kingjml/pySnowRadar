import pytest
from pyWavelet.calcpulse import calcpulsewidth

def test_calcpulsewidth_invalid():
    with pytest.raises(AssertionError):
        assert 2 == 1