import time
#TODO: add look up table for all years
def gpstoutc(gps_time):
    gps_time_ls = _leapsec(gps_time)
    utc_time = [time.gmtime(ds) for ds in gps_time_ls]
    return utc_time

def gpsptoutc(gps_time):
    '''
    This is for the 'packed' style of time that comes from the ATM product
    '''
    #Srem = rem(gpsTimeHHMMSS,1e3)
    #MS = MSrem

    #SSrem = rem(gpsTimeHHMMSS,1e5)
    #SS = (SSrem - MSrem)/1e3

    #MMrem= rem(gpsTimeHHMMSS,1e7)
    #MM = (MMrem - SSrem)/1e5

    #HHrem = rem(gpsTimeHHMMSS,1e9)
    #HH = (HHrem - MMrem)/1e7
    #utcTime = ((HH*3600) + (MM*60) + SS +(MS/1000) - gpsAheadUtc)';
    
    
    
def _leapsec(gps_time):
    '''Correct for leap seconds'''
    leap_epoc = 19
    leap_cur = 36 #2016 leap from http://hpiers.obspm.fr/eop-pc/index.php?index=TAI-UTC_tab&lang=en
    leap_sec = leap_cur - leap_epoc
    gps_time_adj = gps_time - leap_sec
    return gps_time_adj