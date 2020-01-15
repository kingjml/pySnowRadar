import numpy as np
import pywt

def Wavelet_TN(data, null_2_space, delta_fast_time_range, n_snow, ref_snow_layer, cwt_precision):
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
    
    lin_coefs = cwt(data, pywt.Wavelet('haar'), lin_scale_vect, cwt_precision)
    log_coefs = cwt(10 * np.log10(data), pywt.Wavelet('haar'), log_scale_vect, cwt_precision)
    
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
    print(f'locs_as: {locs_as} {type(locs_as)}')
    return locs_as, locs_si

def cwt(data, wavelet, scales, precision):
    '''
    Implementation of the Continuous Wavelet Transform

    Arguments:
        data: preprocessed snowradar signal data
        wavelet: the specific Wavelet to use (currently the Haar wavelet)
        scales: 
        precision: precision to apply to wavelet operations (default 10)

    Outputs:
        out_coefs:(?)
    '''
    out_coefs = np.zeros((np.size(scales), data.size))
    int_psi, x = pywt.integrate_wavelet(wavelet, precision=precision)
    step = x[1] - x[0]
    x_step = (x[-1] - x[0]) + 1
    
    j_a = [np.arange(scale * x_step) / (scale * step) for scale in scales]
    j_m = [np.delete(j, np.where((j >= np.size(int_psi)))[0]) 
           for j in j_a if np.max(j) >= np.size(int_psi)]
    coef_a = [-np.sqrt(scales[i]) 
               * np.diff(np.convolve(data, int_psi[x.astype(np.int)][::-1]))
              for (i, x) in enumerate(j_m)]
    d_a = [(coef.size - data.size) / 2 for coef in coef_a]
    out_coefs = np.asarray([coef[int(np.floor((coef.size - data.size) / 2))
                            :int(-np.ceil((coef.size - data.size) /2))] for coef in coef_a])
    return out_coefs
