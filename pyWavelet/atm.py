import h5py
import numpy as np
import timefunc
from datetime import datetime, timedelta

_qc_pitch_max = 5 #Max ATM pitch in def
_qc_roll_max = 5 #Max ATM roll in deg

class atm:
    def __init__(self, file_path, l_case = 'meta'):
        with h5py.File(file_path, 'r') as f:
            if (l_case=='full'):
                self.elv = f['elevation'][()]
                self.lat = f['latitude'][()]
                self.lon = f['longitude'][()]
                self.packed_time = np.around(f['instrument_parameters']['time_hhmmss'][()],3) #This is a weird packed time format...
                self.elasped = f['instrument_parameters']['rel_time'][()]
                self.pitch = f['instrument_parameters']['pitch'][()]
                self.roll =f['instrument_parameters']['roll'][()]
            self.extent = np.hstack((f['ancillary_data']['min_longitude'][()]-180,
                                     f['ancillary_data']['min_latitude'][()],
                                     f['ancillary_data']['max_longitude'][()]-180,
                                     f['ancillary_data']['max_latitude'][()])).ravel()
                
        self.file_name = file_path.split('\\')[-1]
        self.file_epoch = datetime.strptime(self.file_name [8:23], "%Y%m%d_%H%M%S") 
        self.load_type = l_case
        
        
    def calcutc(self):
        hh_rem = np.mod(self.packed_time, 1e4)
        hh = (self.packed_time-hh_rem)/1e4
        mm = np.modf(hh_rem/1e2)[1]
        ss = np.modf(np.mod(hh_rem, 1e2))[1]
        ms_rem = np.around(np.modf(hh_rem)[0], 3)
        gps_seconds = (hh*3600)+(mm*60)+ss+ms_rem
        gps_ts = self.file_epoch.replace(hour=0, minute=0, second=0, microsecond=0).timestamp()+gps_seconds
        self.utc_time = timefunc.utcleap(gps_ts)