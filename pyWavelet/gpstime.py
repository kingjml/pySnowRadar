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
    
    
def _leapsec(gps_time):
    '''Correct for leap seconds'''
    leap_epoc = 19
    leap_cur = 36 #2016 leap from http://hpiers.obspm.fr/eop-pc/index.php?index=TAI-UTC_tab&lang=en
    leap_sec = leap_cur - leap_epoc
    gps_time_adj = gps_time - leap_sec
    return gps_time_adj