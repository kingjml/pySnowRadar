import numpy as np
RADAR_ROLL_QC = 5 #degrees
RADAR_PITCH_QC = 5 #degrees
N_NOISE_BINS = 100 #Number of noise range bins to sample

def radarqc(data):
    roll_flag = data['Roll'] > RADAR_ROLL_QC
    pitch_flag = data['Roll'] > RADAR_PITCH_QC
    
    #np.apply_along_axis(picklayers.picklayers, 0, data['Data'],n2n,delta_fast_time_range,n_snow)
    
def noisepower(data, noise_stop):
    #noise_power = 
    print(data['Data'])
    #noise_start  = (noise_stop-N_NOISE_BINS)
    #np.apply_along_axis(_noisepowercol, 0, data['Data'], noise_start,noise_stop)
    
    #return np.apply_along_axis(np.median,0,data['Data'][noise_start:noise_stop,:])
    #return np.apply_along_axis(np.median,0,data['Data'])
    
def _noisepowercol(data, n_start, n_stop):
    print(n_start)
    
def atmqc(data):
    print('todo')
    
def snowdepthqc(data):
    print('todo')