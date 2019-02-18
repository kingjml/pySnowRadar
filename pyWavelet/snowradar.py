import os
import h5py
from . import matfunc
from . import timefunc
import numpy as np
from scipy import signal
from shapely.geometry import box

C = 299792458 # Vacuum speed of light
QC_PITCH_MAX = 5 # Max ATM pitch in def
QC_ROLL_MAX = 5 # Max ATM roll in deg

# https://ops.cresis.ku.edu/wiki/index.php/Raw_File_Guide
CRESIS_RAW_FILE_LUT = {
    'snow': 'OIB_MAT',
    'snow2': 'OIB_MAT',
    'snow3': 'OIB_MAT',
    'snow4': 'OIB_MAT',
    'snow5': 'AWI_MAT',
    'snow8': 'OIB_MAT',
    'snow9': 'OIB_MAT',
    'snow10': 'OIB_MAT'
}

# TODO overload the class so it can except botht the .mat and NC snow radar files
class SnowRadar:
    def __init__(self, file_path, l_case):
        self.file_path = os.path.abspath(file_path)
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(self.file_path)
        else:
            print('Processing: ' + self.file_path)
        self.file_name = os.path.basename(self.file_path)
        self.load_type = l_case
        self.air_snow = None
        self.snow_ice = None
        self.epw = None # equiv_pulse_width 
        self.n2n = None # Null to Null space
        # we try to populate these in _populate_instanceattr()
        self.surface = None
        self.elv_corr = None
        # becomes True if FMCW echograms are compressed from source
        self.compressed = False
        radar_dat = matfunc.unified_loader(self)
        self._populate_instanceattr(radar_dat)
        self.epw = None # Equiv pulse width 
        self.n2n = None # Null to Null space
        
    def get_surface(self, smooth=True):
        '''
        Simple surface tracker based on maximum
        This should be refined and is largely a place holder 
        
        '''
        surf_bin = np.argmax(self.data_radar, axis=0)
        surf_time = np.interp(surf_bin,np.arange(0,len(self.time_fast)),self.time_fast)
        return surf_time, surf_bin        
        
    def as_dict(self):
        '''generic method, to be extended by OIB and AWI subclasses'''
        return {
            'fname': self.file_name,
            'fpath': self.file_path,
            'l_case': self.load_type
        }

    def _populate_instanceattr(self, radar_dat):
        '''
        This is part of the unified loader mechanism that automatically accounts
        for the differences between SnowRadar source datasets
        
        Currently-supported source datasets:
        Operation IceBridge         (2016): Matlab v5 HDF
        Operation IceBridge         (2017): Matlab v7 HDF
        Alfred Wegener Institute    (2017): Matlab v7 HDF

        '''
        # scrape some metadata in order to decide how to treat the sourcefile
        radar_name_code = radar_dat['param_records']['radar_name']
        self.data_type = CRESIS_RAW_FILE_LUT.get(radar_name_code, 'Unknown')
        season = radar_dat['param_records']['season_name']
        mission = radar_dat['param_records']['cmd']['mission_names']
        day, segment = radar_dat['param_records']['day_seg'].split('_')
        gps_src = radar_dat['param_records']['gps_source']

        f0 = radar_dat['param_records']['radar']['wfs']['f0']
        f1 = radar_dat['param_records']['radar']['wfs']['f1']
        fmult = radar_dat['param_records']['radar']['wfs']['fmult']
        if self.data_type == 'AWI_MAT':
            # AWI mat files store 2 values for f0, f1 and fmult
            # Here we force the script to use the 2nd value
            f0 = float(f0[1])
            f1 = float(f1[1])
            fmult = float(fmult[1])
        lats = radar_dat['Latitude']
        lons = radar_dat['Longitude']
        gps_times = radar_dat['GPS_time']
        utc_times = [timefunc.utcleap(gps) for gps in gps_times]
        self.file_epoch = utc_times
        fast_times = radar_dat['Time']
        self.bandwidth = np.abs((f1 - f0) * fmult)
        self.dft = fast_times[1] - fast_times[0] # delta fast time
        self.dfr = self.dft * 0.5 * C # delta fast time range
        # load just the metadata concerning timing
        if self.load_type == 'meta':
            time_start = gps_times.min()
            time_end = gps_times.max()
            self.time_gps = np.asarray((time_start, time_end))
            self.time_utc = np.asarray((
                timefunc.utcleap(time_start),
                timefunc.utcleap(time_end)
            ))
        # load the full dataset including as much as possible
        elif self.load_type == 'full':
            data_radar = radar_dat['Data']
            elevation = radar_dat['Elevation']
            if data_radar.shape[0] == elevation.shape[0]:
                self.data_radar = np.transpose(data_radar)
            else:
                self.data_radar = data_radar
            self.elevation = elevation
            self.time_gps = gps_times
            self.time_utc = np.asarray([
                timefunc.utcleap(t) for t in self.time_gps
            ])
            self.time_fast = fast_times
            self.lat = lats
            self.lon = lons
            self.roll = radar_dat['Roll']
            self.pitch = radar_dat['Pitch']     
            # Sometimes the surface is recorded in the matfile
            try:
                self.surface = radar_dat['Surface']
            except KeyError:
                print('Surface unavailable')
            # Sometimes there are elev corrections available in
            # the matfile
            try:
                self.elv_corr = radar_dat['Elevation_Correction'].astype(np.int64)
            except KeyError:
                pass
            # Check if the FMCW echograms are compressed in the matfile
            try:
                self.trunc_bins = radar_dat['Truncate_Bins'].astype(np.int64)
                self.compressed = True
            except KeyError:
                pass
            
            # Check if the file has previously been elevation corrected
            # TODO: check type of output
            try:
                self.elev_corrected = radar_dat['param_records']['get_heights']['elev_correction']
            except KeyError:
                pass
            
        # Geospatial boundary box
        self.extent = np.hstack((
            lons.min(), lats.min(),
            lons.max(), lats.max()
        )).ravel()
        self.poly = box(*self.extent)


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
        # JK: Need to be careful here, f becomes an array if bandwidth is as well.
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
    
    def decompress_data(self):
        '''
        Ported by Josh King from CRESIS uncompress_echogram.m by John Paden
        '''
        print('Decompressing data...') 
        Nz = self.elv_corr.max()
        Nt = Nz + len(self.trunc_bins)
        self.elevation = self.elevation - self.elv_corr * self.dfr
        self.surface = self.surface - self.elv_corr * self.dft
        self.data_radar = np.pad(self.data_radar, [(Nz,0 ), (0, 0)], 'constant', constant_values=(np.nan))
         
        for rline in np.arange(0,self.data_radar.shape[1]):
            self.data_radar[:, rline] = np.roll(self.data_radar[:, rline], -self.elv_corr[rline])
        if len(self.time_fast) == len(self.trunc_bins):
            t0 = self.time_fast[0] - Nz*self.dft
        else:
            t0 =  self.time_fast[self.trunc_bins[0]] - Nz * self.dft
            self.time_fast = t0 + self.dft*np.arange(0, Nt-1)
        print('finished')
    
    def elevation_compensation(self, perm_ice=3.15):
        '''
        Ported by Josh King from CRESIS elevation_compensation.m by John Paden
        https://github.com/kingjml/pyWavelet/blob/master/pyWavelet/legacy/elevation_compensation.m

        Inputs:
            perm_ice: The permittivity of ice (default 3.15)

        Outputs:
            elev_axis: elevation axis based on bin timing and an assumption of permittivity

        Notes:
            This will not work on instances of basic SnowRadar class 
            since they do not have self.data_radar, self.time_utc, 
            self.elevation, self.surface attributes

            MikeB: might be worth investigating use of AbstractBaseClasses/Methods
            or factory methods (https://stackoverflow.com/a/682545)
            for the different inputs (AWI, OIB, CRESIS, NSIDC, etc.)
        '''        
        # Quick check for existance of non-null 'surface' data attribute 
        if self.surface is None:
            print('Elevation compensation not possible without surface estimate')
            return

        # Mikeb: placeholders for commonly-repeated operations
        half_speed_of_light = C * 0.5 
        half_c_in_ice = half_speed_of_light / np.sqrt(perm_ice) 
        time_fast_size = len(self.time_fast)        # TODO: better name for this?

        max_elev = self.elevation.max()
        min_elev = np.min(
            self.elevation - self.surface * half_speed_of_light -
            (self.time_fast[-1] - self.surface) * half_c_in_ice
        )

        # Create an elevation axis based on the bin timing and an assumption of permittivity
        delta_range = self.dft * half_c_in_ice   
        dt_air = delta_range / half_speed_of_light  # TODO: better name for this?
        dt_ice = self.dft                           # TODO: Is this correct?
        elev_axis = np.arange(max_elev, min_elev, -delta_range)

        # Zero-pad the radar data to provide space for interpolation
        zero_pad_len = len(elev_axis) - time_fast_size - 1
        # This will be our new compensated radar data array
        radar_comp = np.concatenate((
            self.data_radar, 
            #np.full((zero_pad_len, self.data_radar.shape[1]), np.nan)),
            np.zeros((zero_pad_len, self.data_radar.shape[1]))),
            axis=0
        )

        # Create the corrections to be applied
        d_range = max_elev - self.elevation
        d_time = d_range / half_speed_of_light
        d_bins = np.round(d_time / dt_ice)

        def create_compensation(row_idx):
            ''' 
            Bad function name because I don't know what it's doing
            
            Can probably be refactored
            
            JK: This function finds the ice interface and scales the time array depending on medium
            '''
            e_val = self.elevation[row_idx]
            s_val = self.surface[row_idx]
            surf_elev = e_val - s_val * half_speed_of_light
            time0 = -(max_elev - e_val) / half_speed_of_light
            last_air_idx = np.max(np.argwhere(elev_axis > surf_elev))
            new_time = (time0 + dt_air * np.arange(0, last_air_idx - 1))
            if last_air_idx < elev_axis.shape[0]:
                first_ice_idx = last_air_idx + 1
                time0 = s_val + (surf_elev - elev_axis[first_ice_idx]) / half_c_in_ice
                new_time = np.concatenate((
                    new_time, 
                    time0 + dt_ice * (
                        np.arange(0, len(elev_axis) - len(new_time) - 1)
                    )), 
                    axis=0
                )
            return new_time

        for idx in range(self.data_radar.shape[1]):
            comp = create_compensation(idx)
            data_subset = self.data_radar[:time_fast_size, idx]
            # ensure that the shapes of self.time_fast and data_subset are the same
            if self.time_fast.shape != data_subset.shape:
                raise BaseException(
                    'Cannot complete elevation compensation: ' +\
                    '\n\ttime_fast and data_subset not same shape at index %s' % idx
                )
            radar_comp[:, idx] = np.interp(
                comp,
                self.time_fast,
                data_subset
            )
        
        # MikeB: I am strongly against having methods that
        # modify class attributes instead of returning outputs.
        #
        # Users will not have any indication that their data
        # has been modified unless they are paying attention or 
        # have read the source code.
        # 
        # There is also no way to backtrack once you run the function,
        # short of re-reading the raw datafiles again.
        self.elevation += d_range
        self.surface += d_time
        self.data_radar = radar_comp
        self.surface_elev = self.elevation - self.surface * half_speed_of_light
        return elev_axis


# The OIB snow radar (2-8 GHz) data comes as matlab v5 or v7
# We use a bit of hacky magic to fit it into numpy arrays from dicts
# TODO: Check file type where NSIDC = NC and CRESIS = MAT
class OIB(SnowRadar):
    def __init__(self, file_path, l_case='meta'):
        super().__init__(file_path, l_case)
        self.data_type = 'OIB_MAT'
        radar_dat = matfunc.unified_loader(self)
        super()._populate_instanceattr(radar_dat)

    def as_dict(self):
        '''extend superclass as_dict method to include more OIB-relevant information'''
        result = super().as_dict()
        result.update({
            'tstart': self.time_utc.min(),
            'tend': self.time_utc.max(),
            'poly': self.poly.wkt,
        })
        return result

    def __str__(self):
        return f'{self.data_type} Datafile: {self.file_name}'


# The AWI snow radar data comes as matlab v7 so its closer to a HDF file
# Use h5py to read and process it
class AWI(SnowRadar):
    def __init__(self, file_path, l_case='meta'):
        super().__init__(file_path, l_case)
        self.data_type = 'AWI_MAT'
        radar_dat = matfunc.unified_loader(self)
        super()._populate_instanceattr(radar_dat)
        
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
        return f'{self.data_type} Datafile: {self.file_name}'