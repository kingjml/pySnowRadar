from pathlib import Path
import numpy as np
import pandas as pd

# Let's shoehorn pathlib in, because I like it better than os
ATM_LOCAL_DIR = Path(__file__).parent / 'data' / 'atm'

QC_PITCH_MAX = 5    # Max allowable pitch (deg)
QC_ROLL_MAX = 5     # Max allowable roll (deg)
QC_DEPTH_MAX = 1.5  # Max allowable snow depth (m) (?)
QC_ELEV_MAX = 0.55  # Max allowable elevation for htopo data (m)
QC_MAX_DIST = 500   # Max allowable distance from nadir (m). For cropping ATM data (?)


def error_check(sr, **kwargs):
    '''
    SnowRadar data integrity check which can be applied before picker processing

    Inputs:
        sr:    A SnowRadar object instance or path to file
    Outputs:
        error_flags:    numpy array containing error flagging indexed as follows
            [0]:  Bad radar amplitude data (0 or NaN value)
            [1]:  Aircraft pitch greater than tollerance
            [2]:  Aircraft roll greater than tollerance
            [3]:  Reserved for future
            [4]:  Reserved for future
            [5]:  Reserved for future
    '''
        
    filled = np.full(sr.data_radar.shape[1], False)
    error_flags = np.stack((~(np.sum(sr.data_radar, axis = 0) > 0), # Check for valid data
                      abs(sr.pitch) > QC_PITCH_MAX, # Check for valid pitch
                      abs(sr.roll) > QC_ROLL_MAX, # Check for valid roll
                      filled, # Reserved flags
                      filled, # Reserved flags
                      filled), # Reserved flags
                      axis = 1)
    return error_flags

def qc_htopo(sr, htopo, snow_depth):
    '''
    ***heavily WIP state***

    Apply quality-control through use of pre-downloaded ATM data and by the following QC flags:
        QC1:    lack of ATM data for given trace
        QC2:    ATM data is available, but too far away from SnowRadar data (?)
        QC3:    pitch or roll outside allowable range
        QC4:    ATM data outside allowable range
        QC5:    Picked Air-Snow interface is "deeper" than Snow-Ice interface
        QC6:    Derived snow depth is outside allowable range

    Inputs:
        sr:             A SnowRadar object instance
        snow_depth:     Derived snow depths from sr
    
    Outputs:
        QC'ed snow_depth

    Based on work by Marissa Dattler (https://www.github.com/mdattler)
    and original MatLab code by Thomas Newman
    '''
    # first thing to check is whether there is any local ATM data on the same day as the sr data
    if not ATM_LOCAL_DIR.exists():
        print(f'***Error: Local ATM directory does not exist: {ATM_LOCAL_DIR.resolve()}')
        return
    atm_data = sorted(ATM_LOCAL_DIR.glob(f'*ATM*_{sr.day}_*.h5'))
    if len(atm_data == 0):
        print(f'***Error: No ATM data found for date: {sr.day}')
        return
    QC3 = np.logical_or(abs(sr.pitch) < QC_PITCH_MAX, abs(sr.roll) < QC_ROLL_MAX)
    