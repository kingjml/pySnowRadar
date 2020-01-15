import numpy as np
from pySnowRadar.algorithms import GSFC_NK, NSIDC, Wavelet_TN

def pick_layers(data, params, picker_func='Wavelet-TN'):
    '''
    General layer-picking function that allows for multiple methods of interface layer detection

    Arguments:
        data: 1D numpy array representing a single snowradar trace
        params: dictionary containing enough key-value pairs for the chosen picker to operate properly
        picker_func: string denoting which picker to use. Currently supporting:
            Wavelet-TN -> default
            NSIDC -> placeholder
            GSFC-NK -> placeholder

    Outputs:
        airsnow_layer: Array index location of where the picker thinks the Air-Snow 
                       interface is for the given snowradar trace
        snowice_layer: Array index location of where the picker thinks the Snow-Ice
                       interface is for the given snowradar trace
    '''
    choice = picker_func.lower()
    if choice == 'wavelet-tn':
        try:            
            airsnow_layer, snowice_layer = Wavelet_TN(
                data, 
                params['n2n'],
                params['dfr'],
                params['n_snow'],
                params.get('ref_snow_layer', 1), # if no ref_snow_layer found, use default of 1
                params.get('cwt_precision', 10) # if no cwt_precision found, use default of 10
            )
        except KeyError:
            raise Exception('Missing or invalid parameters supplied for Wavelet-TN picker')
    elif choice == 'nsidc':
        try:
            airsnow_layer, snowice_layer = NSIDC(
                data,
            )
        except KeyError:
            raise Exception('Missing or invalid parameters supplied for NSIDC picker')
    elif choice == 'gsfc-nk':
        try:
            airsnow_layer, snowice_layer = GSFC_NK(
                data,
            )
        except KeyError:
            raise Exception('Missing or invalid parameters supplied for GSFC-NK picker')
    else:
        raise ValueError('Chosen picker unsupported. Must be one of ["Wavelet_TN", "NSIDC", "GSFC_NK]')
    
    return airsnow_layer, snowice_layer

if __name__ == '__main__':
    print(pick_layers(np.linspace(0,25,25), {'n2n':0}, 'NSIDC'))
