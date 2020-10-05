import logging
from pathlib import Path
import numpy as np
import pandas as pd

# TODO: Need to have a process to init these params
QC_PITCH_MAX = 5    # Max allowable pitch (deg)
QC_ROLL_MAX = 5     # Max allowable roll (deg)
QC_MIN_SNR = 5      # Min SNR at peak (-)

def calc_snr(sr, noise_bins = 100, interface = None):
    '''
    Calculate the signal-to-noise ratio for an entire frame or at 
    specific interfaces
    
    If no interface is provided SNR is calculated as peak-to-noise
    for each trace

    Input:
        noise_bins: The number of leading bins to consider as noise
        interface: Optional numpy array of interface bins to return SNR
    
    Output:
        SNR for an entire frame or at specific bins (interface)
    '''
    noise = np.median(sr.data_radar[0: noise_bins, :], axis = 0)
    if interface is not None:
        with np.errstate(divide='ignore', invalid='ignore'):
            return 10*np.log10(np.array([sr.data_radar[idx, trace] for trace,idx in enumerate(interface)]) / noise)
    else:
        with np.errstate(divide='ignore', invalid='ignore'):
            return 10*np.log10(sr.data_radar.max(axis = 0)/noise)

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
    snr = calc_snr(sr)
    error_flags = np.stack((~(np.sum(sr.data_radar, axis = 0) > 0), # Check for valid data
                      abs(sr.pitch) > QC_PITCH_MAX, # Check for valid pitch
                      abs(sr.roll) > QC_ROLL_MAX, # Check for valid roll
                      np.nan_to_num(snr) < QC_MIN_SNR, # Check for valid SNR
                      filled, # Reserved flags
                      filled), # Reserved flags
                      axis = 1)
    return error_flags