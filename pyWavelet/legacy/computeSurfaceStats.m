function dts = computeSurfaceStats(X,sipId)
% computeSurfaceStats- This function computes the statisitcal
% properties of the input vector X
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    X   - The input vector
%    sipId - calculate surface inflection point flag ('sipYes'/'sipNo')
%
% Outputs:
%    dts - A structure containing the descriptive statisitcs.
%
% Example: 
%    -
%
% Other m-files required:
% - /harmmean.m
% - /kurtosis.m
% - /mad.m
% - /nanmax.m
% - /prctile.m
% - /skewness.m
% - /stats/stats/trimmean.m
% - /wavelet/wavelet/cwt.m
%
% Subfunctions: none
%
% MAT-files required: none
%
% See also: calcSurfaceAttributes
%
%
% Author: Thomas Newman, Ph.D., Earth Sciences
%
% February 2015; Last revision: Base. 

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% BEGIN CODE
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%==========================================================================
%% Measures of Central Tendency.
%==========================================================================
dts.ArithmeticMean = mean(X);
% Arithmetic mean

dts.HarmonicMean = harmmean(X);
% Harmonic mean

dts.TrimmedMean10 = trimmean(X,10); % 10 percent.
% Mean excluding outliers

%==========================================================================
%% Measures of Dispersion.
%==========================================================================
dts.MedianAbsoluteDeviation = mad(X,1);
% The median absolute deviation

dts.MeanAbsoluteDeviation = mad(X,0);
% The mean absolute deviation

dts.StandardDeviation = std(X);
% Standard deviation

%==========================================================================
%% Measures of Shape.
%==========================================================================
dts.Kurtosis = kurtosis(X);
% Kurtosis

dts.Skewness = skewness(X);
% Skewness

p = [0, 1, 2.5, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 97.5, 99, 100];
% Defining the percentiles at which to sample the surface distribution. 

pSummary = prctile(X,p);
% The percentiles used to sample the distribution.

%--------------------------------------------------------------------------
dts.Percentile0000 = pSummary(1);
dts.Percentile0010 = pSummary(2);
dts.Percentile0025 = pSummary(3);
dts.Percentile0050 = pSummary(4);
dts.Percentile0100 = pSummary(5);
dts.Percentile0200 = pSummary(6);
dts.Percentile0300 = pSummary(7);
dts.Percentile0400 = pSummary(8);
dts.Percentile0500 = pSummary(9);
dts.Percentile0600 = pSummary(10);
dts.Percentile0700 = pSummary(11);
dts.Percentile0800 = pSummary(12);
dts.Percentile0900 = pSummary(13);
dts.Percentile0950 = pSummary(14);
dts.Percentile0975 = pSummary(15);
dts.Percentile0990 = pSummary(16);
dts.Percentile1000 = pSummary(17);
%--------------------------------------------------------------------------

%==========================================================================
%% Calculate the surface inflection point.
%==========================================================================
switch sipId
    case 'sipYes'
        
        maxPercentageScale = 10;
        % The maximum wavelet scale defined by a percentage.
        
        maxWaveletScale = length(X)/maxPercentageScale;
        % The maximum  wavelet scale.
        
        scaleVectPre = 2:1:maxWaveletScale;
        % Vector of scales up to the maximum wavelet scale.
        
        scaleVect = scaleVectPre(mod(scaleVectPre,2)==1);
        % Scale vector containing only odd scales.
        
        sortedX = sort(X,'ascend');
        % Sorts the values in input vector X in ascending order.
        
        cwtCoefs = cwt(sortedX,scaleVect,'db1');
        % Computes the continuous wavelet coefficients of the signal vector X
        % at real, positive scales in scaleVect, using wavelet 'haar'.
        
        cwtCoefs(:,1:ceil(maxWaveletScale/2))           = nan;
        cwtCoefs(:,(end-ceil(maxWaveletScale/2)+1:end)) = nan;
        % To negate edge effects the edges are set to nans.
        
        sumCwtCoefs = sum(cwtCoefs,1)./size(cwtCoefs,1);
        % Sums the cwt coefficients and divides by the number of coefs.
        
        [maxSumCwtCoefs,maxSumCwtCoefsInd] = nanmax(sumCwtCoefs);
        % Finds the location of the minimum gradient of the summed CwtCoefs
        
        inflectionPoint = sortedX(maxSumCwtCoefsInd);
        % The value associated with the surface inflection point.
        
        dts.DistibutionInflectionPoint = inflectionPoint;
        % assigning the inflection point to the dts structure.  
end
%==========================================================================

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% END OF CODE
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end