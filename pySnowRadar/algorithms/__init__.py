from .Wavelet import Wavelet_TN, Wavelet_JK
from .GSFC import GSFC_NK, NSIDC

def available_pickers():
    return [Wavelet_TN, Wavelet_JK]