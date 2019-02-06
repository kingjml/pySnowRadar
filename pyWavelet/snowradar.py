import os
import h5py
from . import matfunc
from . import timefunc
import numpy as np
from scipy import signal
from shapely.geometry import box

C = 299792458 #Vacuum speed of light
QC_PITCH_MAX = 5 #Max ATM pitch in def
QC_ROLL_MAX = 5 #Max ATM roll in deg

# TODO overload the class so it can except botht the .mat and NC snow radar files
class SnowRadar:
    def __init__(self, file_path, l_case):
        
        self.file_path = os.path.abspath(file_path)
        self.file_name = os.path.basename(self.file_path)
        self.load_type = l_case
        
        self.air_snow = None
        self.snow_ice = None
        self.epw = None #equiv_pulse_width 
        self.n2n = None #Null to Null space
        
    def as_dict(self):
        '''generic method, to be extended by OIB and AWI subclasses'''
        return {
            'fname': self.file_name,
            'fpath': self.file_path,
            'l_case': self.load_type
        }

    def calcpulsewidth(self, oversample_num=1000, num_nyquist_ts=100):
        '''
        bandwidth: radar bandwidth in hz
        win: window function to be applied to the freq domain
        os_num: the amount to oversample the nyquist by
        
        '''
        # Time Vector
        nyquist_sf = 2 * self.bandwidth
        fs = nyquist_sf * oversample_num 
        time_step = 1 / fs 
        max_time = num_nyquist_ts * oversample_num * time_step
        time_vect = np.arange(-max_time, max_time, time_step)
    
        # Frequency domain object
        half_bandwidth = self.bandwidth / 2
        n_FFT = len(time_vect)
        f = fs * np.linspace(-0.5, 0.5, n_FFT)
        n_band_points = np.sum(np.abs(f) <= half_bandwidth)
    
        # Create spectral window 
        # TODO: Check CRESIS windowing, add options if necessary
        spectral_win = signal.hann(n_band_points)
    
        # Frequency domain processing
        #JK: Need to be careful here, f becomes an array if bandwidth is as well.
        # Change it to use f.shape?
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
        self.epw = equiv_pulse_width_time * C
    
        # Calc null-to-null pulse width
        with np.errstate(divide = 'ignore'):
            invert_l10_power = -10 * np.log10(power_signal_norm)
        peak_idx, _ = signal.find_peaks(invert_l10_power)
        closest_peaks = np.sort(np.abs(peak_idx - max_idx))
    
        null_2_width = 2 * np.mean(closest_peaks[0:1])
        null_2_time = null_2_width * time_step
        self.n2n = null_2_time * C


# The OIB snow radar (2-8 GHz) data comes as matlab v5
# We use a bit of hacky magic to fit it into numpy arrays from dicts
# TODO: Check file type where NSIDC = NC and CRESIS = MAT
class OIB(SnowRadar):
    def __init__(self, file_path, l_case='meta'):
        super().__init__(file_path, l_case)
        radar_dat = matfunc.loadmat(file_path)
        self.data_type = 'OIB_MAT'
        self.bandwidth = np.abs((radar_dat['param_records']['radar']['wfs']['f1'] -
                          radar_dat['param_records']['radar']['wfs']['f0']) *
                          radar_dat['param_records']['radar']['wfs']['fmult'])
        self.dft = radar_dat['Time'][1] - radar_dat['Time'][0]  #delta fast time
        self.dfr = (self.dft / 2) * C #delta fast time range

        lats = radar_dat['Latitude']
        lons = radar_dat['Longitude']
        
        if l_case == 'meta':
            time_start = radar_dat['GPS_time'].min()
            time_end = radar_dat['GPS_time'].max()
            self.time_gps = np.asarray((time_start, time_end))
            self.time_utc = np.asarray((
                timefunc.utcleap(time_start),
                timefunc.utcleap(time_end)
            ))
        elif l_case == 'full':
            self.time_gps = radar_dat['GPS_time']
            self.time_utc = np.asarray([
                timefunc.utcleap(t) for t in self.time_gps
            ])
            self.data_radar = radar_dat['Data']
            self.time_fast = radar_dat['Time']
            self.lat = lats
            self.lon = lons
            self.elevation = radar_dat['Elevation']
            self.roll = radar_dat['Roll']
            self.pitch = radar_dat['Pitch']
            
            if 'Elevation_Correction' in radar_dat.keys():
                self.surface = radar_dat['Surface']
            else:
                print('Surface unavailable')
                self.surface = None
                
            # Looks like there are correcltions available for some files
            # If this field is available its an easy fix
            if 'Elevation_Correction' in radar_dat.keys():
                self.elv_corr = radar_dat['Elevation_Correction']
            else:
                self.elv_corr = None
            
            #Check if the FMCW echograms are compressed at source
            if 'Truncate_Bins' in radar_dat.keys():
                self.compressed = True
                self.trunc_bins = radar_dat['Truncate_Bins']
            else:
                self.compressed = False

        gps_times = radar_dat['GPS_time']
        utc_times = [timefunc.utcleap(gps) for gps in gps_times]
        self.file_epoch = utc_times
        
        self.extent = np.hstack((
            lons.min(), lats.min(),
            lons.max(), lats.max()
        )).ravel()      

        self.poly = box(*self.extent)

    def as_dict(self):
        '''extend superclass as_dict method to include more OIB-relevant information'''
        result = super().as_dict()
        result.update({
            'tstart': self.time_utc.min(),
            'tend': self.time_utc.max(),
            'poly': self.poly.wkt
        })
        return result

    def __str__(self):
        return f'OIB Datafile: {self.file_name}'   


#The AWI snow radar data comes as matlab v7 so its closer to a HDF file
#Use h5py to read and process it
class AWI(SnowRadar):
    def __init__(self, file_path, l_case='meta'):
        super().__init__(file_path, l_case)
        radar_dat = h5py.File(file_path)
        
        self.data_type = 'AWI_MAT'

        #Not sure why but there are two records for f0 and f1. The first one is only a 2Ghz range?
        f0 = radar_dat[list(radar_dat['param_records']['radar']['wfs']['f0'])[1][0]].value
        f1 = radar_dat[list(radar_dat['param_records']['radar']['wfs']['f1'])[1][0]].value
        fmult = radar_dat[list(radar_dat['param_records']['radar']['wfs']['fmult'])[1][0]].value
        self.bandwidth = np.float64(np.abs((f1 - f0) * fmult)) #TODO: The scalar cast is a hack
        
        self.dft = radar_dat['Time'].value[0][1] - radar_dat['Time'].value[0][0]  #delta fast time
        self.dfr = (self.dft / 2) * C #delta fast time range

        lats = radar_dat['Latitude'].value
        lons = radar_dat['Longitude'].value

        if l_case == 'meta':
            time_start = radar_dat['GPS_time'].value.min()
            time_end = radar_dat['GPS_time'].value.max()
            self.time_gps = np.asarray((time_start, time_end))
            self.time_utc = np.asarray((
                timefunc.utcleap(time_start),
                timefunc.utcleap(time_end)
            ))
        
        elif l_case=='full':
            self.time_gps = radar_dat['GPS_time'].value
            self.time_utc = np.asarray([
                timefunc.utcleap(t) for t in self.time_gps])
            self.data_radar = np.transpose(radar_dat['Data'].value)
            self.time_fast = radar_dat['Time'].value.flatten()
            self.lat = lats
            self.lon = lons
            self.elevation = radar_dat['Elevation'].value
            
            # TODO: If there is not a rough surface available
            # generate one quickly for the elevation correction
            if 'Surface' in [key for key in radar_dat.keys()]:
                self.surface = radar_dat['Surface'].value
            else:
                print('Surface unavailable')
                self.surface = None

        self.extent = np.hstack((
            lons.min(), lats.min(),
            lons.max(), lats.max()
        )).ravel()  
        self.poly = box(*self.extent)
        
    def as_dict(self):
        '''extend superclass as_dict method to include more AWI-relevant information'''
        result = super().as_dict()
        result.update({
            'tstart': self.time_utc.min(),
            'tend': self.time_utc.max(),
            'poly': self.poly.wkt
        })
        return result

    def __str__(self):
        return f'AWI Datafile: {self.file_name}'