import h5py
import numpy as np
import datetime

_qc_pitch_max = 5 #Max ATM pitch in def
_qc_roll_max = 5 #Max ATM roll in deg

class atm:
    def __init__(self, file_path, l_case = 'min'):
        with h5py.File(file_path, 'r') as f:
            if (l_case=='full'):
                self.elv = f['elevation'][()]
                self.lat = f['latitude'][()]
                self.lon = f['longitude'][()]
                self.time = np.around(f['instrument_parameters']['time_hhmmss'][()],3)
                self.elasped = f['instrument_parameters']['rel_time'][()]
                self.pitch = f['instrument_parameters']['pitch'][()]
                self.roll =f['instrument_parameters']['roll'][()]
                self.extent = np.hstack((f['ancillary_data']['min_latitude'][()],
                                    f['ancillary_data']['max_latitude'][()],
                                    f['ancillary_data']['min_longitude'][()],
                                    f['ancillary_data']['min_longitude'][()])).ravel()
        self.file_name = file_path.split('\\')[-1]
        self.file_epoch = datetime.datetime.strptime(self.file_name [8:23], "%Y%m%d_%H%M%S") 
        self.load_type = l_case
        
    def calcutc(self):
        hh_rem = np.mod(self.time, 1e4)
        hh = (self.time-hh_rem)/1e4
        mm = np.modf(hh_rem/1e2)[1]
        ss = np.modf(np.mod(hh_rem, 1e2))[1]
        ms_rem = np.around(np.modf(hh_rem)[0], 3)
        gps_seconds = (hh*3600)+(mm*60)+ss+ms_rem
        gps_ts = self.file_epoch.replace(hour=0, minute=0, second=0, microsecond=0).timestamp() + gps_seconds
        self.utc_time = [datetime.datetime.fromtimestamp(ts) for ts in gps_ts]
        
    def _leapsec(self, gps_time):
        '''Correct for leap seconds'''
        leap_epoc = 19
        leap_cur = 36 #2016 leap from http://hpiers.obspm.fr/eop-pc/index.php?index=TAI-UTC_tab&lang=en
        leap_sec = leap_cur - leap_epoc
        gps_time_adj = gps_time - leap_sec
        return gps_time_adj