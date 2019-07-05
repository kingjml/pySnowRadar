import os
import h5py
from datetime import datetime

from . import timefunc

class ATM:
    '''
    Python representation of a NASA Airborne Topographic Mapper (ATM) Laser Altimeter dataset

    Based on work by Marissa Dattler (https://www.github.com/mdattler)
    
    '''
    def __init__(self, file_path):
        self.file_path = os.path.abspath(file_path)
        self.file_name = os.path.basename(self.file_path)
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(self.file_path)
        file_date = self.file_name.split('_')[1]
        self.day = datetime.strptime(file_date, '%Y%m%d')
        ## may be more useful to people
        #file_datetime = ''.join(self.file_name.split('_')[1:3]).split('.')[0]
        #self.start_time = datetime.strptime(file_datetime, '%Y%m%d%H%M%S')
        with h5py.File(self.file_path, mode='r') as src:
            self.data_type = f"ATM {src.attrs['short_name'].decode()}"
            # I hate the [()] notation, but .values is deprecated...
            #     explanation: "[]"" implies an array slice, but since we 
            #     pass an empty slice "()", it returns the full array
            self.pitch = src['instrument_parameters/pitch'][()]
            self.roll = src['instrument_parameters/roll'][()]
            self.latitude = src['latitude'][()]
            # convert longitude from 0:360 to -180:180
            self.longitude = ((src['longitude'][()] - 180) % 360) - 180
            self.elevation = src['elevation'][()]
            # convert time in HHMMSS.ssssss to seconds-since-start-of-day
            self.time_gps = timefunc.atm_hhmmss_to_sec(
               src['instrument_parameters/time_hhmmss'][()]
            )
    
    def __str__(self):
        return f'{self.data_type} Datafile: {self.file_name}'
    