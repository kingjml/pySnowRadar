function [equivalentPulseWidth,null2NullPulseWidth] = calculatePulseWidth(bandwidth,win)
% calculatePulseWidth - This function calculates the equivalent
% pulsewidth and null to null pulsewidth, in space, given an input bandwidth
% and window function
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    bandwidth - The bandwidth of the radar system.
%    win       - The window function applied to the frequency domain.
%
% Outputs:
%    equivalentPulseWidth - This is the equivalent width of the pulse in
%                           space. 
%    null2NullPulseWidth  - This is the null to null pulse width in space.
%
% Example: 
% -
%
% Other m-files required:
% - /signal/signal/hann.m
% - /signal/signal/rectwin.m
% 
%
% Subfunctions: none
%
% MAT-files required: none
%
% See also: pickSnowLayer, qcUncertCalc
%
%
% Author: Thomas Newman, Ph.D., Earth Sciences
%
% February 2015; Last revision: Base. 

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% BEGIN CODE
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%==========================================================================
%% Constants and variables.
%==========================================================================
c = 299792458; 
% The vacuum speed of light.

oversampleNumber = 1000;
% The amount to oversample the nyquist by.

numNyquistTimeSteps = 100;
% The number of nyquist timesteps.

%==========================================================================
%% Constructing the time vector.
%==========================================================================
nyquistSamplingFrequency = 2*bandwidth;
% The nyquist sampling frequency.

fs = nyquistSamplingFrequency*oversampleNumber;
% Sampling frequency.

timeStep = 1/fs;
% The sample spacing in time.

maxTime  = numNyquistTimeSteps*oversampleNumber*timeStep;
% The maximum and minimum time limits of the time vector.

timeVect = -maxTime:timeStep:maxTime;
% The time vector.

%==========================================================================
%% Constructing the frequency domain object.
%==========================================================================
halfBandwidth = bandwidth/2;
% Half the bandwidth.

NFFT          = length(timeVect);
% The number of points to include in the fft.

f             = fs*linspace(-0.5,0.5,NFFT);
% The frequency vector.

numBandPoints = sum(abs(f)<=halfBandwidth);
% The number of points within the bandwidth. 

%==========================================================================
%% Constructing the spectral window.
%==========================================================================
switch win
    case 'rect'
        spectralWin   = rectwin(numBandPoints);
    case 'hann'
        spectralWin   = hann(numBandPoints);
end

%==========================================================================
%% Frequency domain processing.
%==========================================================================
freqDomainSignal                        = zeros(1,length(f));
% Create a vector of zeros to hold the frequency domain signal.

freqDomainSignal(abs(f)<=halfBandwidth) = spectralWin';
% Place the spectral window signal within the frequency domain vector.

shiftFreqDomainSignal                   = ifftshift(freqDomainSignal);
% Inverse FFT shift: swaps the left and right halves of the vector 

timeDomainSignal                        = ifft(shiftFreqDomainSignal)*NFFT;
% Converts frequency domain signal into the time domain.

timeSig                                 = fftshift(timeDomainSignal);
% Shift zero-frequency component to center of spectrum. Moves the 
% zero-frequency component to the center of the array.

powerSignal                             = abs(timeSig.^2);
% Calculates the signal power. 

powerSignalNorm                         = powerSignal./max(powerSignal);
% Calculates the normalised signal power. 

%==========================================================================
%% Calculate the equivalent pulse width.
%==========================================================================
equivalentPulseWidthVal  = sum(powerSignalNorm);
% The equivalent pulse width.

equivalentPulseWidthTime = equivalentPulseWidthVal*timeStep;
% The equivalent pulse width in time.

equivalentPulseWidth     = equivalentPulseWidthTime*c;
% The equivalent pulse width in space.

%==========================================================================
%% Calculate the null-to-null pulse width.
%==========================================================================
[maxVal,maxInd] = max(powerSignalNorm);
% finds the centre of the main lobe.

invertedLog10Power = -10*log10(powerSignalNorm);
% The negative logarithm of the normalised signal power.

[pks,locs]          = findpeaks(invertedLog10Power);
% Find the peaks of the negative logarithm of the normalised signal power.

[closestPeaks,IX]   = sort(abs(locs-maxInd));
% Find the closest peaks to the location of the main lobe.

null2nullWidth      = 2*mean([closestPeaks(1),closestPeaks(2)]);
% Finds the null-to-null width by using the average between the two closest
% peaks in the negative logarithm of the normalised signal power.

null2NullPulseWidthTime  = null2nullWidth*timeStep;
% The null-to-null pulsewidth in time.

null2NullPulseWidth = null2NullPulseWidthTime*c;
% The null-to-null pulsewidth in space.
%==========================================================================

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% END OF CODE
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end