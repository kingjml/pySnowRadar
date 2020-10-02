import numpy as np
from scipy import signal


def Peakiness(data,
              delta_fast_time_range, 
              n_snow,
              log_peak_threshold = 0.7,
              lin_peak_threshold = 0.5, 
              pp_r_threshold = 30, 
              pp_l_threshold = 30,
              **kwargs):
    '''
    Author: Arttu Jutila (AWI)
    
    Function to detect 2 interface layers from a given Snow Radar signal:
        Air-Snow interface
        Snow-Ice Interface

    Arguments:
        data: waveform, 1D radar data array
        delta_fast_time_range: one range bin in meters
        n_snow: the refractive index of snow
        log_peak_threshold: power threshold for a-s interface, default 0.7
        lin_peak_threshold: height threshold for detecting normalized linear peaks, default 0.5
        pp_r_threshold: threshold to filter right-hand peakiness values, default 30
        pp_l_threshold: threshold to filter left-hand peakiness values, default 30

    Outputs:
        locs_as: Air-snow interface bin index (float)
        locs_si: Snow-ice interface bin index (float)

    '''
    try:
        data = data[~np.isnan(data)]
        #normalize data
        data_norm = data / np.nanmax(data)
        #log-scale normalized data
        data_norm_dB = 10 * np.log10(data_norm)
    
        #Calculate noise
        noise_mean = np.nanmean(data_norm_dB[0:100])
    
        #Find logarithmic data peaks that are above the mean noise level by log_peak_threshold
        logpeaks,_ = signal.find_peaks(data_norm_dB, height = noise_mean + log_peak_threshold * (np.nanmax(data_norm_dB) - noise_mean))
    
        #Calculate left-hand peakiness value for the logarithmic data peaks from the linear normalized data
        range_bins = 10 #equals to approximately 2 x 3-dB range resolution (~8 cm)
        PP_L_log = []
        for peak in logpeaks:
            PP_L_log.append(data_norm[peak] / np.mean(data_norm[peak - range_bins:peak]) * range_bins)
    
        #Find where the left-hand peakiness values fulfill the threshold value
        valid_PP_L = np.argwhere(np.array(PP_L_log) >= pp_l_threshold).flatten()
    
        #Find linear data peaks that are above the lin_peak_threshold
        linpeaks,_ = signal.find_peaks(data_norm, height = lin_peak_threshold)
    
        #If more than five peaks in the linear data were detected, the waveform is determined as ambiguous and no interfaces are returned.
        if len(linpeaks) > 5:
            locs_as = np.nan
            locs_si = np.nan
        else:    
            #Calculate right-hand peakiness value for the linear data peaks from the linear normalized data
            PP_R_lin = []
            for peak in linpeaks:
                PP_R_lin.append(data_norm[peak] / np.mean(data_norm[peak + 1:peak + 1 + range_bins]) * range_bins)
    
             #Find where the right-hand peakiness values are valid: filter according to pp_r_threshold, but including the maximum return, and including peaks only within 1.5 m from the estimated air-snow interface
            valid_PP_R = np.argwhere(((np.array(PP_R_lin) >= pp_r_threshold) | 
                                     (linpeaks == np.argmax(data_norm))) & 
                                     (linpeaks < logpeaks[0] + np.ceil(1.5 / (delta_fast_time_range / n_snow)).astype(int))).flatten()
        
            #If no valid peaks after filtering, set interface locations to NaN
            if (len(valid_PP_L) == 0) or (len(valid_PP_R) == 0):
                locs_as = np.nan
                locs_si = np.nan
            else: #Filter peak locations
                logpeaks = logpeaks[valid_PP_L]
                linpeaks = linpeaks[valid_PP_R]
                                 
                #Make sure that the assumed air-snow interface comes before the snow-ice interface
                if (logpeaks[0] > linpeaks[-1]):
                    locs_as = np.nan
                    locs_si = np.nan
                else:
                    #Air-snow interface is located at the first valid logarithmic data peak
                    locs_as = 1.0 * logpeaks[0]
                    #Snow-ice interface is located at the last valid linear data peak
                    locs_si = 1.0 * linpeaks[-1]
    except:
        locs_as = np.nan
        locs_si = np.nan
    
    return locs_as, locs_si