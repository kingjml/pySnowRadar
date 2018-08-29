import h5py
import numpy as np
import timefunc
from datetime import datetime, timedelta, timezone
from shapely.geometry import box

_qc_pitch_max = 5 #Max ATM pitch in def
_qc_roll_max = 5 #Max ATM roll in deg

class atm:
    def __init__(self, file_path, l_case = 'meta'):
        with h5py.File(file_path, 'r') as f:
            self.time_packed = np.around(f['instrument_parameters']['time_hhmmss'][()],3)[[0,-1]] #This is a weird time format...
            self.extent = np.hstack((np.mod((f['ancillary_data']['min_longitude'][()]+180),360)-180,
                                     f['ancillary_data']['min_latitude'][()],
                                     np.mod((f['ancillary_data']['max_longitude'][()]+180),360)-180,
                                     f['ancillary_data']['max_latitude'][()])).ravel()
            if (l_case=='full'):
                self.elv = f['elevation'][()]
                self.lat = f['latitude'][()]
                self.lon = np.mod((f['longitude'][()]+180),360)-180
                #self.elasped = f['instrument_parameters']['rel_time'][()]
                self.pitch = f['instrument_parameters']['pitch'][()]
                self.roll =f['instrument_parameters']['roll'][()]
                self.time_packed = np.around(f['instrument_parameters']['time_hhmmss'][()],3) #This is a weird time format...
            
        self.poly = box(self.extent[0], self.extent[1], self.extent[2], self.extent[3])        
        self.file_name = file_path.split('\\')[-1]
        self.file_path = file_path
        self.load_type = l_case
        
        #Sort out the file times
        self.time_epoch = datetime.strptime(self.file_name [8:23], "%Y%m%d_%H%M%S").replace(tzinfo=timezone.utc).timestamp()
        self.calcutc()
        
        
    def as_dict(self):
        if (self.load_type=='meta'):
            return {'fname': self.file_name, 
                    'fpath': self.file_path, 
                    'tstart':self.time_utc[0], 
                    'tend':self.time_utc[-1],
                    'poly': self.poly}
        
    def calcutc(self):
        hh_rem = np.mod(self.time_packed, 1e4)
        hh = (self.time_packed-hh_rem)/1e4
        mm = np.modf(hh_rem/1e2)[1]
        ss = np.modf(np.mod(hh_rem, 1e2))[1]
        ms_rem = np.around(np.modf(hh_rem)[0], 3)
        gps_seconds = (hh*3600)+(mm*60)+ss+ms_rem
        time_start = datetime.utcfromtimestamp(self.time_epoch)
        gps_ts = time_start.replace(hour=0, minute=0, second=0, microsecond=0, tzinfo=timezone.utc).timestamp()+gps_seconds
        self.time_utc = timefunc.utcleap(gps_ts)