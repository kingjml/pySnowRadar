from .Wavelet import Wavelet_TN, Wavelet_JK
from .GSFC import GSFC_NK, NSIDC # Not yet integrated
from .Peakiness import Peakiness

def available_pickers():
    return [Wavelet_TN, Wavelet_JK, Peakiness]