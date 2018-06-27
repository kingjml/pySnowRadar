function snow = pickSnowLayer(frame)
% pickSnowLayer - This function picks the snow layer of the current
% frame using wavelet techniques and outputs a structure of snow depths and
% ancillary information.
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    frame - This is a structure containing the currently loaded frame.
%
% Outputs:
%    snow - This is a structure containing the picked snow depths and
%    ancillary information.
%
% Example:
% -
%
% Other m-files required:
% - /wavelet/wavelet/cwt.m
% - calculatePulseWidth
% 
% Subfunctions: none
%
% MAT-files required: none
%
% See also: waveletSnowProcessor
%
%
% Author: Thomas Newman, Ph.D., Earth Sciences
%
% February 2015; Last revision: Base.

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% BEGIN CODE
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%==========================================================================
% Material properties.
%==========================================================================
densitySnow = 0.32;
% Alexandrov, V. et al., 2010
% The relation between sea ice thickness and freeboard in the Arctic
% The Cryosphere, 4(3), pp.373-380.

permSnow = (1+0.51*densitySnow)^3;
% Panzer, B. et al., 2013.
% An ultra-wideband, microwave radar for measuring snow thickness on
% sea ice and mapping near-surface internal layers in polar firn.
% Journal of Glaciology, 59(214), pp.244-254.

n_snow = sqrt(permSnow);
% To cover from permitivity to refractive index.

c = 299792458;
% The vacuum speed of light.

%==========================================================================
% Radar properties.
%==========================================================================
bandwidth               = frame.radarBandwidth;
% The exact bandwidth of the radar system.

deltaFastTime           = frame.Time(2) - frame.Time(1);
% The delta fast time.

deltaFastTimeRange      = (deltaFastTime/2) *c;
% The delta fast time converted to vaccuum range.

%==========================================================================
% Wavelet scale properties.
%==========================================================================

%--------------------------------------------------------------------------
% To pick the snow/ice interface - Linear space.
%--------------------------------------------------------------------------
[pulseWidth,null2nullPulseWidth] =  calculatePulseWidth(bandwidth,'hann');
% calculates the pulsewidth and null2null pulsewidth

refScaleLinMetres = 2*null2nullPulseWidth;
% -1 -1 -1 -1 -1  0  1  1  1  1  1
%                   [null2nullPulseWidth]
% The pulse width should be entirely contained to the right of the zeros.

maxScaleLin = ceil(refScaleLinMetres/deltaFastTimeRange);
% The maximum wavelet scale to detect the snow/ice interface.

linearScaleVectPre = 2:1:maxScaleLin;
% All scales.

linearScaleVect = linearScaleVectPre(mod(linearScaleVectPre,2)==1);
% Only odd scales.

%--------------------------------------------------------------------------
% To pick the air/snow interface - Logarithmic space.
%--------------------------------------------------------------------------
referenceSnowLayer = 1;
% 100 cm reference snow layer.
% Constrained by time constraints.

snowLayerOpl         = referenceSnowLayer*n_snow*2;
% The optical path length of two-way travel through the reference snow layer.

refScaleLogMetres    = 2*snowLayerOpl;
% -1 -1 -1 -1 -1  0  1  1  1  1  1
%                    [snowLayerOpl]
% The referenceSnowLayer width should be entirely contained to the right of the zeros.

maxScaleLog = ceil(refScaleLogMetres/deltaFastTimeRange);
% The maximum wavelet scale to detect the air/snow interface.

logarithmicScaleVectPre = 2:1:maxScaleLog;
% All scales.

logarithmicScaleVect = logarithmicScaleVectPre(mod(logarithmicScaleVectPre,2)==1);
% Only odd scales.

%--------------------------------------------------------------------------
% Preallocation to speed up processing.
%--------------------------------------------------------------------------
switch frame.radarType
    
    case 'Snow_Radar_Deconvolved'
        rangeBinVectAirSnow = zeros(1,frame.echogramNumberOfCols);
        % The air/snow preallocation vector.
        
        rangeBinVectSnowIce = zeros(1,frame.echogramNumberOfCols);
        % The snow/ice preallocation vector.
        
    case 'Snow_Radar'
        rangeBinVectAirSnow = nan(1,frame.echogramNumberOfCols);
        % The air/snow preallocation vector.
        
        rangeBinVectSnowIce = nan(1,frame.echogramNumberOfCols);
        % The snow/ice preallocation vector.
        
end
%--------------------------------------------------------------------------

%==========================================================================
% The Wavelet snow depth picking algorithm.
%==========================================================================
switch frame.radarType
    
    case 'Snow_Radar_Deconvolved'
        
        for iTrace = 1:frame.echogramNumberOfCols;
            
            %----------------------------------------------------------------------
            % Snow/ice range bin - Linear waveform - Picked first.
            %----------------------------------------------------------------------
            dataCol = double(frame.Data(:,iTrace));
            % The current trace converted from single to double precision.
            
            linearCoefs = cwt(dataCol,linearScaleVect,'db1');
            % Computes the continuous wavelet coefficients of the signal vector dataCol
            % at real, positive scales, using wavelet 'haar'.
            
            linearCoefs(:,1:ceil(maxScaleLin/2))           = 0;
            linearCoefs(:,(end-ceil(maxScaleLin/2)+1:end)) = 0;
            % To negate edge effects.
            
            sumLinearCoefs = sum(linearCoefs,1)./size(linearCoefs,1);
            % Sums the linear coefficients and divides by size.
            
            [pksSi,locsSi] = max(-sumLinearCoefs);
            % Finds the maximum of the summed linear coefficients.
            
            if ~isempty(pksSi)
                snowIceRangeBin = locsSi;
                % The snow/ice range bin is the maximum of the summed linear
                % coefficients,
                
                rangeBinVectSnowIce(iTrace)  = snowIceRangeBin;
                % Places the snowIceRangeBin for output.
                
            else
                rangeBinVectSnowIce(iTrace)  = nan;
                % If no maximum sumLinearCoefs then places a nan for output.
                
            end
            
            %----------------------------------------------------------------------
            % Air/snow range bin - Log waveform - Picked second
            %----------------------------------------------------------------------
            dataColLog = 10*log10(dataCol);
            % The current trace in logarithmic space.
            
            logarithmicCoefs = cwt(dataColLog,logarithmicScaleVect,'db1');
            % Computes the continuous wavelet coefficients of the signal vector dataCol
            % at real, positive scales, using wavelet 'haar'.
            
            logarithmicCoefs(:,1:ceil(maxScaleLog))           = 0;
            logarithmicCoefs(:,(end-ceil(maxScaleLog)+1:end)) = 0;
            % To negate edge effects.
            
            sumLogarithmicCoefs = sum(logarithmicCoefs,1)./size(logarithmicCoefs,1);
            % Sums the logarithmic coefficients and divides by size.
            
            [pksAs,locsAs] = max(-sumLogarithmicCoefs);
            % Finds the maximum of the summed logarithmic coefficients.
            
            if ~isempty(pksAs)
                
                airSnowRangeBin = locsAs;
                % The air/snow range bin is the maximum of the summed
                % logarithmic coefficients.
                
                rangeBinVectAirSnow(iTrace)  = airSnowRangeBin;
                % Places the airSnowRangeBin in a vector for output.
                
            else
                rangeBinVectAirSnow(iTrace)  = nan;
                % If no maximum sumLogarithmicCoefs then places a nan for output.
            end
            %----------------------------------------------------------------------
            
        end
end

%==========================================================================
% Converting the range bin separation into a snow depth.
%==========================================================================
snowDepth = ((rangeBinVectSnowIce-rangeBinVectAirSnow)*deltaFastTimeRange)/n_snow;
% The derived snow depth.

%==========================================================================
% Outputs.
%==========================================================================

%--------------------------------------------------------------------------
% Input properties.
%--------------------------------------------------------------------------
snow.frameLatitude  = frame.Latitude;
snow.frameLongitude = frame.Longitude;
snow.frameGpsTime   = frame.GPS_time;
snow.frameFastTime  = frame.Time;

snow.radarType      = frame.radarType;
snow.dateStr        = frame.dateStr;
snow.segmentNumber  = frame.segmentNumber;
snow.frameNumber    = frame.frameNumber;
snow.fullFilePath   = frame.fullFilePath;

%--------------------------------------------------------------------------
% Time properties.
%--------------------------------------------------------------------------
snow.frameUtcSecs       = frame.utcSecs;
snow.frameMinUtcSecs    = frame.minUtcSecs;
snow.frameMaxUtcSecs    = frame.maxUtcSecs;
snow.frameMedianUtcSecs = frame.medianUtcSecs;

%--------------------------------------------------------------------------
% Echogram properties.
%--------------------------------------------------------------------------
snow.echogramNumberOfRows = frame.echogramNumberOfRows;
snow.echogramNumberOfCols = frame.echogramNumberOfCols;

%--------------------------------------------------------------------------
% Radar properties.
%--------------------------------------------------------------------------
snow.frameRadarBandwidth = frame.radarBandwidth;

%--------------------------------------------------------------------------
% Material properties.
%--------------------------------------------------------------------------
snow.density             = densitySnow;
snow.permittivity        = permSnow;
snow.refractiveIndex     = n_snow;
snow.c                   = c;

%--------------------------------------------------------------------------
% Wavelet properties.
%--------------------------------------------------------------------------
snow.refScaleLinMetres   = null2nullPulseWidth;
snow.refScaleLogMetres   = refScaleLogMetres;

snow.maxLinearScale      = max(linearScaleVect);
snow.maxLogarithmicScale = max(logarithmicScaleVect);

%--------------------------------------------------------------------------
% Snow layer properties.
%--------------------------------------------------------------------------
snow.rangeBinVectSnowIce = rangeBinVectSnowIce;
snow.rangeBinVectAirSnow = rangeBinVectAirSnow;

snow.depth               = snowDepth;

%--------------------------------------------------------------------------
% The time the file was processed.
%--------------------------------------------------------------------------
snow.dateFrameLoaded     = frame.dateLastLoaded;
snow.dateProcessed       = datestr(clock);
snow.frameTrackingNumber = frame.trackingNumber;

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% END OF CODE
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end