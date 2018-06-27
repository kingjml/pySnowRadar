import numpy as np
import pywt

def picklayers(data,null_2_space,delta_fast_time_range, n_snow, ref_snow_layer = 1):
    ref_scale_lin_m = 2*null_2_space
    max_scale_lin = np.ceil(ref_scale_lin_m/delta_fast_time_range)
    lin_scale_vect_pre = np.arange(2,max_scale_lin,1)
    lin_scale_vect = lin_scale_vect_pre [np.mod(lin_scale_vect_pre,2)==1]
    
    snow_layer_opl = ref_snow_layer*n_snow*2
    ref_scale_log_m = 2*snow_layer_opl
    max_scale_log = np.ceil(ref_scale_log_m/delta_fast_time_range)
    log_scale_vect_pre = np.arange(2,max_scale_log,1)
    log_scale_vect = log_scale_vect_pre [np.mod(log_scale_vect_pre,2)==1]
    
    lin_coefs = cwt(data, pywt.Wavelet('haar'),lin_scale_vect)
    log_coefs = cwt(10*np.log10(data), pywt.Wavelet('haar'),log_scale_vect)
    
    sum_log_coefs = np.sum(log_coefs,axis=0)/log_coefs.shape[0]
    sum_lin_coefs = np.sum(lin_coefs,axis=0)/lin_coefs.shape[0]
    
    locs_si = np.argmax(-sum_lin_coefs)
    locs_as = np.argmax(-sum_log_coefs)

    return locs_as,locs_si

def cwt(data, wavelet, scales, precision = 10):
    out_coefs = np.zeros((np.size(scales), data.size))
    for i in np.arange(np.size(scales)):
        int_psi, x = pywt.integrate_wavelet(wavelet, precision=precision)
        step = x[1] - x[0]
        j = np.floor(
            np.arange(scales[i] * (x[-1] - x[0]) + 1) / (scales[i] * step))
        if np.max(j) >= np.size(int_psi):
            j = np.delete(j, np.where((j >= np.size(int_psi)))[0])
        coef = - np.sqrt(scales[i]) * np.diff(
            np.convolve(data, int_psi[j.astype(np.int)][::-1]))
        d = (coef.size - data.size) / 2.
        out_coefs[i, :] = coef[int(np.floor(d)):int(-np.ceil(d))]
    return out_coefs