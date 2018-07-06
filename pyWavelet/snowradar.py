import matfunc
import timefunc
import numpy as np
from scipy import signal

_c = 299792458
_leap_epoc = 19 #Leap seconds at GPS epoch

#TODO overload the class so it can except botht the .mat and NC snow radar files
class snowradar:
    def __init__(self, file_path, l_case = 'meta'):
        radar_dat =  matfunc.loadmat(file_path)
        self.gps_time = radar_dat['GPS_time']
        self.utc_time = timefunc.utcleap(radar_dat['GPS_time'])
        self.bandwidth = np.abs((radar_dat['param_records']['radar']['wfs']['f1']-
                          radar_dat['param_records']['radar']['wfs']['f0'])*
                          radar_dat['param_records']['radar']['wfs']['fmult'])
        self.dft = radar_dat['Time'][1] - radar_dat['Time'][0]  #delta fast time
        self.dfr = (self.dft/2)*_c #detla fast time range
        self.load_type = l_case
        
        self.air_snow = None
        self.snow_ice = None
        self.epw = None #equiv_pulse_width 
        self.n2n = None #Null to Null space
        
        self.extent = np.hstack((radar_dat['Longitude'].min(),
                                 radar_dat['Latitude'].min(),
                                 radar_dat['Longitude'].max(),
                                 radar_dat['Latitude'].max())).ravel()
                                 
        if (l_case=='full'):
            self.data = radar_dat['Data']
            self.ft = radar_dat['Time']
            self.lat = radar_dat['Latitude']
            self.lon  = radar_dat['Latitude']
        
    def calcpulsewidth(self, window = 'hann', oversample_num = 1000, num_nyquist_ts = 100):
        '''
        bandwidth: radar bandwidth in hz
        win: window function to be applied to the freq domain
        os_num: the amount to oversample the nyquist by
        '''
        # Time Vector
        nyquist_sf = 2*self.bandwidth
        fs = nyquist_sf *oversample_num 
        time_step = 1/fs 
        max_time  = num_nyquist_ts*oversample_num*time_step
        time_vect = np.arange(-max_time,max_time,time_step)
    
        # Frequency domain object
        half_bandwidth = self.bandwidth/2
        n_FFT= len(time_vect)
        f = fs*np.linspace(-0.5,0.5,n_FFT)
        n_band_points = np.sum(np.abs(f)<=half_bandwidth)
    
        # Create spectral window 
        # TODO: allow other windows i.e. rect
        spectral_win = signal.hann(n_band_points)
    
        # Frequency domain processing
        freq_domain_signal = np.zeros(len(f))
        freq_domain_signal[np.abs(f)<half_bandwidth] = spectral_win
        shift_freq_domain_signal = np.fft.ifftshift(freq_domain_signal)
        time_domain_signal = np.fft.ifft(shift_freq_domain_signal)*n_FFT
        time_sig = np.fft.fftshift(time_domain_signal)
        power_signal = np.abs(time_sig**2)
        power_signal_norm = power_signal/np.max(power_signal)
        max_val = np.max(power_signal_norm)
        max_idx = np.argmax(power_signal_norm)
    
        # Calc the equivalent pulse width
        equiv_pulse_width_val = np.sum(power_signal_norm)
        equiv_pulse_width_time = equiv_pulse_width_val*time_step
        self.epw = equiv_pulse_width_time*_c
    
        # Calc null-to-null pulse width
        with np.errstate(divide='ignore'):
            invert_l10_power = -10*np.log10(power_signal_norm)
        peak_idx, _ = signal.find_peaks(invert_l10_power)
        closest_peaks =  np.sort(np.abs(peak_idx - max_idx))
    
        null_2_width = 2*np.mean(closest_peaks[0:1])
        null_2_time = null_2_width*time_step
        self.n2n = null_2_time*_c