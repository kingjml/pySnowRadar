import numpy as np
import pywt

def Wavelet_JK(data, null_2_space, delta_fast_time_range, n_snow, ref_snow_layer, cwt_precision):
    '''
    Function to detect 2 interface layers from a given SnowRadar signal:
        Air-Snow interface
        Snow-Ice Interface

    Currently uses the Continuous Wavelet Transform (cwt) method originally developed
    by Thomas Newman

    Arguments:
        data: 1D radar data array
        null_2_space:(?)
        delta_fast_time_range:(?)
        n_snow: the refractive index of snow(?)
        ref_snow_layer:(?) (default 1)
        cwt_precision: precision arg for cwt (default 10)

    Outputs:
        locs_as: Air-snow interface bin index (integer)
        locs_si: Snow-ice interface bin index (integer)

    '''
    ref_scale_lin_m = 2 * null_2_space
    max_scale_lin = np.ceil(ref_scale_lin_m / delta_fast_time_range)
    lin_scale_vect = np.arange(2, max_scale_lin, 1)[1::2]
    
    snow_layer_opl = ref_snow_layer * n_snow * 2
    ref_scale_log_m = 2 * snow_layer_opl
    max_scale_log = np.ceil(ref_scale_log_m / delta_fast_time_range)
    log_scale_vect = np.arange(2, max_scale_log, 1)[1::2]
    
    lin_coefs, lin_freq = pywt.cwt(data,lin_scale_vect,'gaus1')
    log_coefs, log_freq = pywt.cwt(10 * np.log10(data),lin_scale_vect,'gaus1')

    #lin_coefs = cwt(data, pywt.Wavelet('haar'), lin_scale_vect, cwt_precision)
    #log_coefs = cwt(10 * np.log10(data), pywt.Wavelet('haar'), log_scale_vect, cwt_precision)
    
    # Negating edge effects here, we use half the max scale on either end
    # Some discussion is needed on this approach because it can sometimes lead to weird picks
    # TODO: Explore other filtering methods (Signal windowing?)
    # TODO: Refactor to support np.nan instead
    lin_coefs[:, 0:np.ceil(max_scale_lin/2).astype(int)] = 0
    lin_coefs[:, -np.ceil(max_scale_lin/2).astype(int):] = 0

    log_coefs[:, 0:np.ceil(max_scale_log/2).astype(int)] = 0
    log_coefs[:, -np.ceil(max_scale_log/2).astype(int):] = 0
    
    sum_log_coefs = np.sum(log_coefs,axis=0) / log_coefs.shape[0]
    sum_lin_coefs = np.sum(lin_coefs,axis=0) / lin_coefs.shape[0]
    
    locs_si = np.argmax(-sum_lin_coefs)
    locs_as = np.argmax(-sum_log_coefs)
    
    return locs_as, locs_si
    
