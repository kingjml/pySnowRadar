function snowQc = qcUncertCalc(frame,surfaceAttributes,snow)
% qcUncertCalc - Performs quality control on the picked snow depths
% using radar system parameters and surface characteristics.
% Assigns a precision for each snow depth pick that passes the quality
% control criterion otherwise assigns a nan.
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    frame             - A structure containing the frame data.
%    surfaceAttributes - A structure containing the surfaceAttributes.
%    snow              - A structure containing the snow depth data.
%
% Outputs:
%    snowQc - A structure containing the quality controlled snow depths
%            and ancillary data.
%
% Example:
% -
%
% Other m-files required:
% - /mad.m
% - calculatePulseWidth
%
% Subfunctions: none
%
% MAT-files required: none
%
% Other files required: 
%
% See also: waveletSnowProcessor
%
%
% Author: Thomas Newman, Ph.D., Earth Sciences
%
% February 2016; Last revision: Base.

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% BEGIN CODE
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%==========================================================================
%% Assign struct to variables.
%==========================================================================
Data                = frame.Data;
% The snow radar echogram.

Time                = frame.Time;
% The frame time vector.


rangeBinVectAirSnow = snow.rangeBinVectAirSnow;
% The air/snow range bin vector.

rangeBinVectSnowIce = snow.rangeBinVectSnowIce;
% The snow/ice range bin vector.

snowDepth           = snow.depth;
% The snow depth vector.


meanAtmPitchVect    = surfaceAttributes.meanAtmPitch;
% The vector of the aircraft pitch, each one calculated from the 200 ATM
% points closest to nadir.

meanAtmRollVect     = surfaceAttributes.meanAtmRoll;
% The vector of the mean ATM rolls, each one calculated from the 200 ATM
% points closest to nadir.

h_topoVect          = surfaceAttributes.h_topo;
% The h_topoVect vector.

%--------------------------------------------------------------------------
% ATM data availability.
%--------------------------------------------------------------------------
atmDataAvailable    = ~isnan(h_topoVect);
% Is ATM data available for the trace.

anyAtmDataAvailable = any(atmDataAvailable);
% Is any ATM data available.

%--------------------------------------------------------------------------
% Snow depth data availability.
%--------------------------------------------------------------------------
snowDataAvailable  = ~isnan(snowDepth);
% Are snow data available for the trace.

anySnowDataAvailable = any(snowDataAvailable);
% Is any snow data available.

%==========================================================================
%% Radar parameters.
%==========================================================================
[pulseWidth,null2nullPulseWidth] = calculatePulseWidth(frame.radarBandwidth,'hann');
% Calculates the pulsewidth and null2null pulsewidth.

rangeResolution = pulseWidth/2;
% The radar range resolution.

rangeResolutionInSnow = rangeResolution/snow.refractiveIndex;
% The range resolution in snow.

deltaFastTime   = Time(2) - Time(1);
% Range bin width: time.

deltaFastTimeRange = (deltaFastTime/2) *snow.c;
% Range bin width: range.

numOfNoiseRangeBins  = 100;
% The number of range bins from which to extract the noise power statisitcs.


interfaceOrderVect = sign(rangeBinVectSnowIce-rangeBinVectAirSnow);
% 1 = snow-ice interface below air snow.

%==========================================================================
%% Quality controls.
%==========================================================================

%--------------------------------------------------------------------------
% QC2.
%--------------------------------------------------------------------------
minDistAtmPoint2NadirQc = 5;
% The distance to the closest ATM point to nadir must be
% <= minDistAtmPoint2NadirQc (5 metres).

%--------------------------------------------------------------------------
% QC3.
%--------------------------------------------------------------------------
maxDistAtmPoint2NadirQc = 50;
% The distance to the 200th closest ATM point to nadir must be
% <= maxDistAtmPoint2NadirQc (50 metres).

%--------------------------------------------------------------------------
% QC4.
%--------------------------------------------------------------------------
meanAtmPitchQc = 5; % degrees.
% The aircraft pitch must be <= meanAtmPitchQc (5 degrees).

%--------------------------------------------------------------------------
% QC5.
%--------------------------------------------------------------------------
meanAtmRollQc  = 5; % degrees.
% The aircraft roll must be <= meanAtmRollQc (5 degrees).

%--------------------------------------------------------------------------
% QC6.
%--------------------------------------------------------------------------
h_topoQc = 0.5; % 0.5
% The h_topoVect value must be <= 0.5 m

%--------------------------------------------------------------------------
% QC7.
%--------------------------------------------------------------------------
interfaceOrderQc = 1;
% The interface order must be 1 as this indicates the snow/ice interface is
% below the air/snow interface.

%--------------------------------------------------------------------------
% QC8.
%--------------------------------------------------------------------------
minSnowDepthQc = rangeResolutionInSnow;
% The minimum detectable snow depth.

%--------------------------------------------------------------------------
% QC9.
%--------------------------------------------------------------------------
maxSnowDepthQc = 2.0; % Changed from 1.5 to 2.0 m
% The snow depth must be <= maxSnowDepth (1.5 metres).

%--------------------------------------------------------------------------
% QC10.
%--------------------------------------------------------------------------
minSnrAirSnowQc = 1;
% The minimum signal-to-noise ratio at air/snow interface must be
% >= minSnrAirSnowQc

%--------------------------------------------------------------------------
% QC11.
%--------------------------------------------------------------------------
minSnrSnowIceQc = 1;
% The minimum signal-to-noise ratio at snow/ice interface must be
% >= minSnrSnowIceQc

%==========================================================================
%% Filter snow depths using quality control criteria and assign uncertainty
%==========================================================================
numQcParameters = 12;
% The number of qc parameters.

qcArray  = ones(numQcParameters,frame.echogramNumberOfCols);
% Prealllocated the Qc array.

noisePowerVect            = nan(1,frame.echogramNumberOfCols);
noisePowerMadVect         = nan(1,frame.echogramNumberOfCols);
noiseAndClutterPowerVect  = nan(1,frame.echogramNumberOfCols);%
signalPowerAirSnowVect    = nan(1,frame.echogramNumberOfCols);
signalPowerSnowIceVect    = nan(1,frame.echogramNumberOfCols);
snrAirSnowVect            = nan(1,frame.echogramNumberOfCols);%
snrSnowIceVect            = nan(1,frame.echogramNumberOfCols);%
noisePowerMadSigmaVect    = nan(1,frame.echogramNumberOfCols);
sixSigmaThresholdVect     = nan(1,frame.echogramNumberOfCols);
noisePowerThresholdVect   = nan(1,frame.echogramNumberOfCols);%
minDistAtmPoint2NadirVect = nan(1,frame.echogramNumberOfCols);% 
maxDistAtmPoint2NadirVect = nan(1,frame.echogramNumberOfCols);%
% Vector preallocation.


%%
%--------------------------------------------------------------------------
% The loop.
%--------------------------------------------------------------------------
for ind = 1:frame.echogramNumberOfCols
    
    %----------------------------------------------------------------------
    % The echogram trace.
    %----------------------------------------------------------------------
    echogramTrace        = Data(:,ind);
    % The current echogram trace.
    
    %----------------------------------------------------------------------
    % Calculating the noise power.
    %----------------------------------------------------------------------
    % For the deconvolved snow radar data the noise power is
    % defined as the average (median) of 100 range bins above
    % the picked air/snow interface.
    
    %......................................................................
    % Corrects for air-snow range bin being to close to top of echogram.
    %......................................................................
    if anySnowDataAvailable==1
        
        if rangeBinVectAirSnow(ind) <= (numOfNoiseRangeBins+1)
            
            noisePowerRangeBinVector = 1:numOfNoiseRangeBins;
            % The range bins containing the noise power.
            
        elseif rangeBinVectAirSnow(ind) > (numOfNoiseRangeBins+1)
            
            noisePowerRangeBinVector = ...
                (rangeBinVectAirSnow(ind)-numOfNoiseRangeBins):1:(rangeBinVectAirSnow(ind)-1);
            % The range bins containing the noise power.
        end
        
    end
    
    %......................................................................
    % Noise power,
    %......................................................................
    if anySnowDataAvailable==1
        
        noisePowerI = median(echogramTrace(noisePowerRangeBinVector));
        % The noise power value.
        
        noisePowerVect(ind) = noisePowerI;
        % The noise power vector.
        
    elseif anySnowDataAvailable==0
        
        noisePowerVect(ind) = nan;
        % The noise power vector.
        
    end
    
    %......................................................................
    % Noise power MAD.
    %......................................................................
    if anySnowDataAvailable == 1
        
        noisePowerMadI = mad(echogramTrace(noisePowerRangeBinVector),1);
        % The noise floor mad value.
        
        noisePowerMadVect(ind) = noisePowerMadI;
        % The noise power mad vector.
        
    elseif anySnowDataAvailable == 0
        
        noisePowerMadVect(ind) = nan;
        % The noise power mad vector.
        
    end
    
    %......................................................................
    % Noise power MAD sigma.
    %......................................................................
    if anySnowDataAvailable == 1
        
        noisePowerMadSigmaI = 1.4826*noisePowerMadI;
        % For normally distributed data, multiply mad by one of the following
        % factors to obtain an estimate of the normal scale parameter sigma.
        
        noisePowerMadSigmaVect(ind) = noisePowerMadSigmaI;
        % The noise power mad sigma vector.
        
    elseif anySnowDataAvailable == 0
        
        noisePowerMadSigmaVect(ind) = nan;
        % The noise power mad sigma vector.
        
    end
    
    %......................................................................
    % Six sigma threshold.
    %......................................................................
    if anySnowDataAvailable == 1
        
        sixSigmaThresholdI = 6*noisePowerMadSigmaI;
        % Standard threshold is 6 sigma.
        
        sixSigmaThresholdVect(ind) = sixSigmaThresholdI;
        % The six sigma threshold vect.
        
    elseif anySnowDataAvailable == 0
        
        sixSigmaThresholdVect(ind) = nan;
        % The six sigma threshold vect.
        
    end
    
    %......................................................................
    % Noise power threshold.
    %......................................................................
    if anySnowDataAvailable == 1
        
        noisePowerThresholdI = noisePowerI+sixSigmaThresholdI;
        % The noise power threshold.
        
        noisePowerThresholdVect(ind) = noisePowerThresholdI;
        % The noise power threshold vect.
        
    elseif anySnowDataAvailable == 0
        
        noisePowerThresholdVect(ind) = nan;
        % The noise power threshold vect.
        
    end
    
    %----------------------------------------------------------------------
    % Calculating the signal power.
    %----------------------------------------------------------------------
    
    %......................................................................
    % Air/snow interface.
    %......................................................................
    if anySnowDataAvailable == 1
        
        signalPowerAirSnowI  = echogramTrace(rangeBinVectAirSnow(ind));
        % The signal power for the air/snow interface pick.
        
        signalPowerAirSnowVect(ind) = signalPowerAirSnowI;
        % The signal power for the air/snow interface pick vector.
        
        snrAirSnowI          = signalPowerAirSnowI/noisePowerI;
        % The SNR for the air/snow interface pick.
        
        snrAirSnowVect(ind)  = snrAirSnowI;
        % The snr at the snow/ice interface pick.
        
    elseif anySnowDataAvailable == 0
        
        snrAirSnowVect(ind)  = nan;
        % The snr at the snow/ice interface pick.
        
    end
    
    %......................................................................
    % Snow/ice interface.
    %......................................................................
    if anySnowDataAvailable == 1
        
        signalPowerSnowIceI  = echogramTrace(rangeBinVectSnowIce(ind));
        % The signal power for the snow/ice interface pick.
        
        signalPowerSnowIceVect(ind) = signalPowerSnowIceI;
        % The signal power for the snow/ice interface pick vector.
        
        snrSnowIceI          = signalPowerSnowIceI/noisePowerI;
        % The SNR for the snow/ice interface pick.
        
        snrSnowIceVect(ind)  = snrSnowIceI;
        % The snr at the snow/ice interface pick.
        
    elseif anySnowDataAvailable == 0
        
        signalPowerSnowIceVect(ind) = nan;
        % The signal power for the snow/ice interface pick vector.
        
        snrSnowIceVect(ind)  = nan;
        % The snr at the snow/ice interface pick.
        
    end
    
    %----------------------------------------------------------------------
    % Calculating the clutter + noise power.
    %----------------------------------------------------------------------
    % The clutter + noise power is defined as the average power in range bins
    % from one range bin after the air/snow pick to
    % one range bin before the snow/ice pick
    %
    % median((a|s)+1:1:(s|i)-1).
    
    if anySnowDataAvailable == 1
        
        noiseAndClutterPowerRangeBinVector = ...
            (rangeBinVectAirSnow(ind)+1):1:(rangeBinVectSnowIce(ind)-1);
        % The range bins containing the noise and clutter power.
        
        noiseAndClutterPowerI = median(echogramTrace(noiseAndClutterPowerRangeBinVector));
        % The noise and clutter power value.
        
        noiseAndClutterPowerVect(ind) = noiseAndClutterPowerI;
        % The noise and clutter power vector.
        
    elseif anySnowDataAvailable == 0
        
        noiseAndClutterPowerVect(ind) = ind;
        % The noise and clutter power vector.
        
    end
    
    %----------------------------------------------------------------------
    % QC1: Quality control by the existence of ATM data.
    %----------------------------------------------------------------------
    if atmDataAvailable(ind) == 0;
        qcArray(1,ind) = 0;
    end
    
    % If a NaN is assigned here, it means that for that particular trace
    % the number of ATM points in a search radius of 500 m about nadir is
    % less than the required number of ATM Points (200).
    
    %----------------------------------------------------------------------
    % QC2: Quality control by minimum distance.
    %----------------------------------------------------------------------
    if atmDataAvailable(ind) == 1;
        
        minDistAtmPoint2Nadir = surfaceAttributes.dtsAtmDist{1, ind}.Percentile0000;
        % The distance to the closest ATM point to nadir.
        
        %..................................................................
        % If fails qc condition.
        %..................................................................
        if minDistAtmPoint2Nadir > minDistAtmPoint2NadirQc
            qcArray(2,ind) = 0;
        end
        %..................................................................
        
        minDistAtmPoint2NadirVect(ind) = minDistAtmPoint2Nadir;
        
    elseif atmDataAvailable(ind) == 0;
        
        minDistAtmPoint2NadirVect(ind) = nan;
        qcArray(2,ind) = 0;
        
    end
    %----------------------------------------------------------------------
    % QC3: Quality control by maximum distance.
    %----------------------------------------------------------------------
    if atmDataAvailable(ind) == 1;
        
        maxDistAtmPoint2Nadir = surfaceAttributes.dtsAtmDist{1, ind}.Percentile1000;
        % The distance to the 200th closest ATM point to nadir.
        
        %..................................................................
        % If fails qc condition.
        %..................................................................
        if maxDistAtmPoint2Nadir > maxDistAtmPoint2NadirQc
            qcArray(3,ind) = 0;
        end
        %..................................................................
        
        maxDistAtmPoint2NadirVect(ind) = maxDistAtmPoint2Nadir;
        
    elseif atmDataAvailable(ind) == 0;
        
        maxDistAtmPoint2NadirVect(ind) = nan;
        qcArray(3,ind) = 0;
        
    end
    %----------------------------------------------------------------------
    % QC4: Quality control by aircraft pitch.
    %----------------------------------------------------------------------
    if atmDataAvailable(ind) == 1;
        
        %..................................................................
        % If fails qc condition.
        %..................................................................
        if abs(meanAtmPitchVect(ind)) > meanAtmPitchQc
            qcArray(4,ind) = 0;
        end
        %..................................................................
        
    elseif atmDataAvailable(ind) == 0;
        
        qcArray(4,ind) = 0;
        
    end
    %----------------------------------------------------------------------
    % QC5: Quality control by aircraft roll.
    %----------------------------------------------------------------------
    if atmDataAvailable(ind) == 1;
        
        %..................................................................
        % If fails qc condition.
        %..................................................................
        if abs(meanAtmRollVect(ind)) > meanAtmRollQc
            qcArray(5,ind) = 0;
        end
        %..................................................................
        
    elseif atmDataAvailable(ind) == 0;
        
        qcArray(5,ind) = 0;
        
    end
    %----------------------------------------------------------------------
    % QC6: Quality control by h_topoVect.
    %----------------------------------------------------------------------
    if atmDataAvailable(ind) == 1;
        
        %..................................................................
        % If fails qc condition.
        %..................................................................
        if h_topoVect(ind) > h_topoQc
            qcArray(6,ind) = 0;
        end
        %..................................................................
        
    elseif atmDataAvailable(ind) == 0;
        
        qcArray(6,ind) = 0;
        
    end
    %----------------------------------------------------------------------
    % QC7: Quality control by Interface order.
    %----------------------------------------------------------------------
    if snowDataAvailable(ind) == 1;
        
        %..................................................................
        % If fails qc condition.
        %..................................................................
        if interfaceOrderVect(ind) < interfaceOrderQc;
            qcArray(7,ind) = 0;
        end
        %..................................................................
        % air/snow range bin must be < snow/ice range bin.
        
    elseif snowDataAvailable(ind) == 0;
        
        qcArray(7,ind) = 0;
        
    end
    %----------------------------------------------------------------------
    % QC8: Quality control by minimum snow depth.
    %----------------------------------------------------------------------
    if snowDataAvailable(ind) == 1;
        
        %..................................................................
        % If fails qc condition.
        %..................................................................
        if snowDepth(ind) < minSnowDepthQc;
            qcArray(8,ind) = 0;
        end
        %..................................................................
        
    elseif snowDataAvailable(ind) == 0;
        
        qcArray(8,ind) = 0;
        
    end
    %----------------------------------------------------------------------
    % QC9: Quality control by maximum snow depth.
    %----------------------------------------------------------------------
    if snowDataAvailable(ind) == 1;
        
        %..................................................................
        % If fails qc condition.
        %..................................................................
        if snowDepth(ind) > maxSnowDepthQc;
            qcArray(9,ind) = 0;
        end
        %..................................................................
        
    elseif snowDataAvailable(ind) == 0;
        
        qcArray(9,ind) = 0;
        
    end
    %----------------------------------------------------------------------
    % QC10: Quality control by signal-to-noise ratio at air/snow.
    %----------------------------------------------------------------------
    if snowDataAvailable(ind) == 1;
        
        %......................................................................
        % If fails qc condition.
        %......................................................................
        if snrAirSnowI < minSnrAirSnowQc
            qcArray(10,ind) = 0;
        end
        %......................................................................
        
    elseif snowDataAvailable(ind) == 0;
        
        qcArray(10,ind) = 0;
        
    end
    %----------------------------------------------------------------------
    % QC11: Quality control by signal-to-noise ratio at snow/ice.
    %----------------------------------------------------------------------
    if snowDataAvailable(ind) == 1;
        
        %......................................................................
        % If fails qc condition.
        %......................................................................
        if snrSnowIceI < minSnrSnowIceQc
            qcArray(11,ind) = 0;
        end
        %......................................................................
        
    elseif snowDataAvailable(ind) == 0;
        
        qcArray(11,ind) = 0;
        
    end
    %----------------------------------------------------------------------
    % QC12: Quality control by clutter noise power.
    %----------------------------------------------------------------------
    if snowDataAvailable(ind) == 1;
        
        %......................................................................
        % If fails qc condition.
        %......................................................................
        if noiseAndClutterPowerVect(ind) < noisePowerThresholdVect(ind)
            qcArray(12,ind) = 0;
        end
        %......................................................................
        
    elseif snowDataAvailable(ind) == 0;
        
        qcArray(12,ind) = 0;
        
    end
    %----------------------------------------------------------------------
    
end

%==========================================================================
%% Calculating the radar range uncertainties.
%==========================================================================
sampUncertAirSnow = deltaFastTimeRange * ones(1,frame.echogramNumberOfCols);
% Sampling uncertainty: air/snow interface.

sampUncertSnowIce = deltaFastTimeRange * ones(1,frame.echogramNumberOfCols);
% Sampling uncertainty: snow/ice interface.

snUncertAirSnow   = rangeResolution./sqrt(2*(signalPowerAirSnowVect./noisePowerVect));
% Signal-to-noise uncertainty: air/snow interface.

sncUncertSnowIce  = rangeResolution./sqrt(2*(signalPowerSnowIceVect./noiseAndClutterPowerVect));
% Signal-to-noise+clutter uncertainty: snow/ice interface.

radarPrecision = sqrt(sampUncertAirSnow.^2 + sampUncertSnowIce.^2 +...
    snUncertAirSnow.^2 + sncUncertSnowIce.^2);
% The calculated snow pick uncertainty.

%==========================================================================
%% Calculating the surface parameters.
%==========================================================================
elevationMinVect  = nan(1,frame.echogramNumberOfCols);
elevation05Vect   = nan(1,frame.echogramNumberOfCols);
elevation50Vect   = nan(1,frame.echogramNumberOfCols);
elevation95Vect   = nan(1,frame.echogramNumberOfCols);
elevationMaxVect  = nan(1,frame.echogramNumberOfCols);
elevationMeanVect = nan(1,frame.echogramNumberOfCols);
elevationStdVect  = nan(1,frame.echogramNumberOfCols);
numAtmPointsVect  = nan(1,frame.echogramNumberOfCols);
surfElevIpVect    = nan(1,frame.echogramNumberOfCols);

for indx = 1:frame.echogramNumberOfCols
    
    if isstruct(surfaceAttributes.dtsAtmElev{1, indx}) == 1
        elevationMinVect(indx)  = surfaceAttributes.dtsAtmElev{1, indx}.Percentile0000;
        elevation05Vect(indx)   = surfaceAttributes.dtsAtmElev{1, indx}.Percentile0050;
        elevation50Vect(indx)   = surfaceAttributes.dtsAtmElev{1, indx}.Percentile0500;
        elevation95Vect(indx)   = surfaceAttributes.dtsAtmElev{1, indx}.Percentile0950;
        elevationMaxVect(indx)  = surfaceAttributes.dtsAtmElev{1, indx}.Percentile1000;
        elevationMeanVect(indx) = surfaceAttributes.dtsAtmElev{1, indx}.ArithmeticMean;
        elevationStdVect(indx)  = surfaceAttributes.dtsAtmElev{1, indx}.StandardDeviation;
        %surfElevIpVect(indx)    = surfaceAttributes.dtsAtmElev{1, indx}.DistibutionInflectionPoint;
        numAtmPointsVect(indx)  = surfaceAttributes.numAtmPointsForStats;
        
    elseif isstruct(surfaceAttributes.dtsAtmElev{1, indx}) == 0
        elevationMinVect(indx)  = nan;
        elevation05Vect(indx)   = nan;
        elevation50Vect(indx)   = nan;
        elevation95Vect(indx)   = nan;
        elevationMaxVect(indx)  = nan;
        elevationMeanVect(indx) = nan;
        elevationStdVect(indx)  = nan;
        %surfElevIpVect(indx)    = nan;
        numAtmPointsVect(indx)  = nan;
    end
    
end

%==========================================================================
%% For intercomparison exercise.
%==========================================================================
timeTruncate = frame.Time(frame.Truncate_Bins);
% The truncated time vector.

deltaRangeBinTime = timeTruncate - timeTruncate(1);
% Delta range bin time.

deltaRangeBinOnewayRange = deltaRangeBinTime*snow.c;
% The delta range bin converted to one-way range.

fastTimeRange     = deltaRangeBinOnewayRange(2) - deltaRangeBinOnewayRange(1);
fastTimeRangeVect = fastTimeRange + zeros(1,length(rangeBinVectAirSnow));
% The fast time range vector.

rangeToAirSnowInterface = zeros(1,length(snowDepth));
rangeToSnowIceInterface = zeros(1,length(snowDepth));
rangeBetweenInterfaces  = zeros(1,length(snowDepth));

for ind2 = 1:length(snowDepth)
    
    if  ~isnan(snowDepth(ind2))
        
        rangeToAirSnowInterface(ind2) = deltaRangeBinOnewayRange(rangeBinVectAirSnow(ind2));
        % The range to the air-snow interface.
        
        rangeToSnowIceInterface(ind2) = deltaRangeBinOnewayRange(rangeBinVectSnowIce(ind2));
        % The range to the snow-ice interface.
        
        rangeBetweenInterfaces(ind2)  = abs(rangeToSnowIceInterface(ind2)-rangeToAirSnowInterface(ind2));
        % The range distance between the two interfaces.
        
    elseif  isnan(snowDepth(ind2))
        
        rangeToAirSnowInterface(ind2) = nan;
        % The range to the air-snow interface.
        
        rangeToSnowIceInterface(ind2) = nan;
        % The range to the snow-ice interface.
        
        rangeBetweenInterfaces(ind2)  = nan;
        % The range distance between the two interfaces.
        
    end
    
    
end

%==========================================================================
%% Output.
%==========================================================================

%--------------------------------------------------------------------------
% From frame structure.
%--------------------------------------------------------------------------
snowQc.frameGpsTime         = frame.GPS_time;
snowQc.frameFastTime        = frame.Time;

snowQc.radarType            = frame.radarType;
snowQc.dateStr              = frame.dateStr;
snowQc.segmentNumber        = frame.segmentNumber;
snowQc.frameNumber          = frame.frameNumber;
snowQc.fullFilePath         = frame.fullFilePath;

snowQc.frameMinUtcSecs      = frame.minUtcSecs;
snowQc.frameMaxUtcSecs      = frame.maxUtcSecs;
snowQc.frameMedianUtcSecs   = frame.medianUtcSecs;

snowQc.dateFrameLoaded      = frame.dateLastLoaded;
snowQc.frameTrackingNumber  = frame.trackingNumber;

snowQc.echogramNumberOfRows = frame.echogramNumberOfRows;
snowQc.echogramNumberOfCols = frame.echogramNumberOfCols;

snowQc.frameRadarBandwidth  = frame.radarBandwidth;

%--------------------------------------------------------------------------
% Book keeping.
%--------------------------------------------------------------------------
snowQc.dateNumVect      = zeros(1,frame.echogramNumberOfCols) + str2num(frame.dateStr);
snowQc.segmentNumVect   = zeros(1,frame.echogramNumberOfCols) + frame.segmentNumber;
snowQc.frameNumVect     = zeros(1,frame.echogramNumberOfCols) + frame.frameNumber;
snowQc.traceNumVect     = 1:frame.echogramNumberOfCols;

%--------------------------------------------------------------------------
% Snow depth.
%--------------------------------------------------------------------------
snowQc.rangeBinVectAirSnow = rangeBinVectAirSnow;
% The air/snow range bin vector.

snowQc.rangeBinVectSnowIce = rangeBinVectSnowIce;
% The snow/ice range bin vector.

snowQc.snowDepth      = snowDepth;
% The snow depth vector.

snowQc.radarPrecision = radarPrecision;
% The snow depth uncertainty.

%--------------------------------------------------------------------------
% Power parameters.
%--------------------------------------------------------------------------
snowQc.noisePowerVect         = noisePowerVect;
snowQc.noisePowerMadVect      = noisePowerMadVect;     

snowQc.signalPowerAirSnowVect = signalPowerAirSnowVect; 
snowQc.signalPowerSnowIceVect = signalPowerSnowIceVect; 

snowQc.noisePowerMadSigmaVect = noisePowerMadSigmaVect; 
snowQc.sixSigmaThresholdVect  = sixSigmaThresholdVect;  

%--------------------------------------------------------------------------
% Quality controls.
%--------------------------------------------------------------------------
snowQc.qcArray = qcArray;
% The qc array.

snowQc.atmDataAvailable  = atmDataAvailable;
% Assign to struct.

snowQc.snowDataAvailable = snowDataAvailable;
% Assign to struct.

%..........................................................................
% QC2
%..........................................................................
snowQc.minDistAtmPoint2NadirVect = minDistAtmPoint2NadirVect;
% A vector of minimum distances to nadir.

snowQc.minDistAtmPoint2NadirQc   = minDistAtmPoint2NadirQc;
% The minimum distance qc conditon.

%..........................................................................
% QC3
%..........................................................................
snowQc.maxDistAtmPoint2NadirVect = maxDistAtmPoint2NadirVect;
% A vector of maximum distances to nadir.

snowQc.maxDistAtmPoint2NadirQc   = maxDistAtmPoint2NadirQc;
% The maximum distance qc conditon.

%..........................................................................
% QC4
%..........................................................................
snowQc.meanAtmPitchVect = meanAtmPitchVect;
% A vector of the mean ATM pitch.

snowQc.meanAtmPitchQc   = meanAtmPitchQc;
% The mean ATM pitch qc conditon.

%..........................................................................
% QC5
%..........................................................................
snowQc.meanAtmRollVect = meanAtmRollVect;
% A vector of the mean ATM roll.

snowQc.meanAtmRollQc   = meanAtmRollQc;
% The mean ATM roll qc conditon.

%..........................................................................
% QC6
%..........................................................................
snowQc.h_topoVect     = h_topoVect;
% The h_topoVect vector.

%..........................................................................
% QC7
%..........................................................................
snowQc.interfaceOrderVect = interfaceOrderVect;
% The interface order vector.

%..........................................................................
% QC8
%..........................................................................
snowQc.minSnowDepthQc = minSnowDepthQc;
% The minimum allowable snow depth.

%..........................................................................
% QC9
%..........................................................................
snowQc.maxSnowDepthQc = maxSnowDepthQc;
% The maximum allowable snow depth.

%..........................................................................
% QC10
%..........................................................................
snowQc.snrAirSnowVect = snrAirSnowVect;
% The SNR at the air/snow interface.

snowQc.minSnrAirSnowQc = minSnrAirSnowQc;
% The minimum condition.

%..........................................................................
% QC11
%..........................................................................
snowQc.snrSnowIceVect = snrSnowIceVect;
% The SNR at the snow/ice interface.

snowQc.minSnrSnowIceQc = minSnrSnowIceQc;
% The minimum condition.

%..........................................................................
% QC12
%..........................................................................
snowQc.noiseAndClutterPowerVect = noiseAndClutterPowerVect;
% The noise and clutter power vector.

snowQc.noisePowerThresholdVect  = noisePowerThresholdVect;
% The noise power threshold vector.

%--------------------------------------------------------------------------
% Surface parameters.
%--------------------------------------------------------------------------
snowQc.elevationMinVect  = elevationMinVect;
snowQc.elevation05Vect   = elevation05Vect;
snowQc.elevation50Vect   = elevation50Vect;
snowQc.elevation95Vect   = elevation95Vect;
snowQc.elevationMaxVect  = elevationMaxVect;
snowQc.elevationMeanVect = elevationMeanVect;
snowQc.elevationStdVect  = elevationStdVect;
snowQc.numAtmPointsVect  = numAtmPointsVect;
snowQc.surfElevIpVect    = surfElevIpVect;

%==========================================================================
% For intercomparison.
%==========================================================================
snowQc.frameUtcSecs            = frame.utcSecs;   % The frame utc seconds.
snowQc.frameLatitude           = frame.Latitude;  % The frame latitude.
snowQc.frameLongitude          = frame.Longitude; % The frame longitude.

snowQc.rangeBetweenInterfaces  = rangeBetweenInterfaces;  % The range between interfaces.
snowQc.rangeToAirSnowInterface = rangeToAirSnowInterface; % The range to the air-snow interface.
snowQc.rangeToSnowIceInterface = rangeToSnowIceInterface; % The range to the snow-ice interface.

snowQc.fastTimeRange           = fastTimeRangeVect;   % The fast time range.
%==========================================================================

snowQc.dateProcessed        = datestr(clock);

%==========================================================================

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% END OF CODE
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end