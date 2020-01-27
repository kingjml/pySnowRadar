import numpy as np

def NSIDC(data):
    '''
    This is a placeholder for the NASA Goddard product hosted at NSIDC.

    Arguments must include the snowradar trace data itself (passed as a 1D float array) as well as 
    any parameters required by the algorithm for layer-picking

    Only 2 picked layers are expected
        1. air-snow layer pick (integer indexing the passed snowradar trace)
        2. snow-ice layer pick (integer indexing the passed snowradar trace)
    
    '''
    as_pick = np.int64(0)
    si_pick = np.int64(0)
    return as_pick, si_pick
