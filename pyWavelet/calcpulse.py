import numpy as np
from scipy import signal
C = 299792458 #Vacuum speed of light

def calcpulsewidth(bandwidth, window='hann', oversample_num=1000, num_nyquist_ts=100):
    '''
    bandwidth: radar bandwidth in hz
    win: window function to be applied to the freq domain
    os_num: the amount to oversample the nyquist by
    '''
    # Time Vector
    nyquist_sf = 2 * bandwidth
    fs = nyquist_sf * oversample_num 
    time_step = 1 / fs 
    max_time  = num_nyquist_ts * oversample_num * time_step
    time_vect = np.arange(-max_time, max_time, time_step)
    
    # Frequency domain object
    half_bandwidth = bandwidth / 2
    n_FFT= len(time_vect)
    f = fs * np.linspace(-0.5, 0.5, n_FFT)
    n_band_points = np.sum(np.abs(f) <= half_bandwidth)
    
    # Create spectral window 
    # TODO: allow other windows i.e. rect
    spectral_win = signal.hann(n_band_points)
    
    # Frequency domain processing
    freq_domain_signal = np.zeros(len(f))
    freq_domain_signal[np.abs(f) < half_bandwidth] = spectral_win
    shift_freq_domain_signal = np.fft.ifftshift(freq_domain_signal)
    time_domain_signal = np.fft.ifft(shift_freq_domain_signal) * n_FFT
    time_sig = np.fft.fftshift(time_domain_signal)
    power_signal = np.abs(time_sig ** 2)
    power_signal_norm = power_signal / np.max(power_signal)
    max_val = np.max(power_signal_norm)
    max_idx = np.argmax(power_signal_norm)
    
    # Calc the equivalent pulse width
    equiv_pulse_width_val = np.sum(power_signal_norm)
    equiv_pulse_width_time = equiv_pulse_width_val * time_step
    equiv_pulse_width = equiv_pulse_width_time * C
    
    # Calc null-to-null pulse width
    with np.errstate(divide='ignore'):
        invert_l10_power = -10 * np.log10(power_signal_norm)
    peak_idx, _ = signal.find_peaks(invert_l10_power)
    closest_peaks = np.sort(np.abs(peak_idx - max_idx))
    
    null_2_width = 2 * np.mean(closest_peaks[0:1])
    null_2_time = null_2_width * time_step
    null_2_space = null_2_time * C
    
    return equiv_pulse_width, null_2_space