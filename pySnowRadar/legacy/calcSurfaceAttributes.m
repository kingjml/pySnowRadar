function surfaceAttributes = calcSurfaceAttributes(frame)
% calcSurfaceAttributes - This function calculates the surface
% attributes using ATM data pertaining to the currently loaded frame.
% These surface attributes are used for quality control and topographic
% filtering, primarily using the h_topo term.
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    frame - This is a structure containing the currently loaded frame.
%
% Outputs:
%    surfaceAttributes - A structure containing properties of the ice
%    surface derived from ATM data.
%
% Example:
% -
%
% Other m-files required:
% - computeSurfaceStats; 
% - getAtmDataForFrame; 
% - latLon2IbXy;
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
%% Setup.
%==========================================================================
numTraces = frame.echogramNumberOfCols;
% The number of traces in the currently loaded snow radar file.

atmTimeCrop = 3; % Seconds (3).
% The number of seconds to crop the ATM data to past the maximum and
% minimum frame utc seconds.

nadirDistanceCrop = 500; % Metres (500).
% The distance from nadir to crop the ATM data to, to increase speed.

numAtmPointsForStats = 200; % Integer (200).
% The number of ATM points from which to calculate surface statistics.

%==========================================================================
%% Load ATM data.
%==========================================================================
atm = getAtmDataForFrame(frame);
% Loads the ATM data applicable for the current snow radar
% frame.

atmFrameCropLower = atm.utcTime >= (frame.minUtcSecs - atmTimeCrop);
atmFrameCropUpper = atm.utcTime <= (frame.maxUtcSecs + atmTimeCrop);
atmFrameCropInd   = atmFrameCropLower & atmFrameCropUpper;
% Creates ATM data mask for data that is ± atmTimeCrop of
% the frame edges.

utcSecsAtmFrameCrop     = atm.utcTime(atmFrameCropInd);
latAtmFrameCrop         = atm.latitude(atmFrameCropInd);
lonAtmFrameCrop         = atm.longitude(atmFrameCropInd);
elevationAtmFrameCrop   = atm.elevation(atmFrameCropInd);
pitchAtmFrameCrop       = atm.pitch(atmFrameCropInd);
rollAtmFrameCrop        = atm.roll(atmFrameCropInd);
% Crops the ATM data so that it is only ± atmTimeCrop of
% the frame edges.

%==========================================================================
%% Convert latitude and longitude to polar stereographic xy.
%==========================================================================

%--------------------------------------------------------------------------
% Converts the frame latitude and longitude to polar stereographic xy coordinates.
%--------------------------------------------------------------------------
[xFrame, yFrame] = latLon2IbXy(frame.Latitude,frame.Longitude);

%--------------------------------------------------------------------------
% Converts the cropped ATM data latitude and longitude to polar stereographic xy coordinates.
%--------------------------------------------------------------------------
if ~isempty(latAtmFrameCrop)
    [xAtmFrameCrop, yAtmFrameCrop] = latLon2IbXy(latAtmFrameCrop,lonAtmFrameCrop);
    % If ATM data exists then perform the conversion.
    
elseif isempty(latAtmFrameCrop)
    xAtmFrameCrop = [];
    yAtmFrameCrop = [];
    % If no ATM data exists then output an empty vector.
end

%==========================================================================
%% Preallocation to increase speed.
%==========================================================================
meanAtmPointUtcSecs = nan(1,numTraces);
meanAtmPointLat     = nan(1,numTraces);
meanAtmPointLon     = nan(1,numTraces);
meanAtmPitch        = nan(1,numTraces);
meanAtmRoll         = nan(1,numTraces);

h_topo              = nan(1,numTraces);
floeElevWav         = nan(1,numTraces);

dtsAtmDist          = cell(1,numTraces);
dtsAtmElev          = cell(1,numTraces);

%==========================================================================
%% The loop.
%==========================================================================
for ind = 1:numTraces;
    % For each trace in the currently loaded snow radar file.
    
    if ~isempty(latAtmFrameCrop)
        % If there is positional data for the current trace.
        
        %------------------------------------------------------------------
        % Filter by distance.
        %------------------------------------------------------------------
        dist = hypot(xFrame(ind)-xAtmFrameCrop,yFrame(ind)-yAtmFrameCrop);
        % Calculates the horizontal distance between the location of the
        % current trace and the cropped ATM points.
        
        distInd = find(dist<nadirDistanceCrop);
        % Keep only the ATM points within nadirDistanceCrop of nadir
        % (nominally 500 m).
        
        [sortDist,sortDistIndPre] = sort(dist(distInd),'ascend');
        % Sort the horizontal distances between the location of the
        % current trace and the cropped ATM points in ascending order.
        
        sortDistInd = distInd(sortDistIndPre);
        % Applies the sorting criteria.
        
        %------------------------------------------------------------------
        % If the number of local ATM points is greater than, or equal to,
        % the required number of ATM Points (nominally 200).
        %------------------------------------------------------------------
        if length(sortDist)>=numAtmPointsForStats
            
            closestAtmPointsInd = sortDistInd(1:numAtmPointsForStats);
            % The indexes of the closest ATM points to nadir (nominally 200).
            
            %..............................................................
            % ATM: UTC times, latitudes and longitudes.
            %..............................................................
            meanAtmPointUtcSecs(ind) = mean(utcSecsAtmFrameCrop(closestAtmPointsInd));
            % The mean UTC times associated with the points closest to
            % nadir (UTC seconds).
            
            meanAtmPointLat(ind)     = mean(latAtmFrameCrop(closestAtmPointsInd));
            % The mean latitudes associated with the points closest to
            % nadir (decimal degrees).
            
            meanAtmPointLon(ind)     = mean(lonAtmFrameCrop(closestAtmPointsInd));
            % The mean longitudes associated with the points closest to
            % nadir (decimal degrees).
            
            %..............................................................
            % ATM: Aircraft attitude.
            %..............................................................
            meanAtmPitch(ind)       = mean(pitchAtmFrameCrop(closestAtmPointsInd));
            % The mean aircraft pitch associated with the points closest to
            % nadir (decimal degrees).
            
            meanAtmRoll(ind)        = mean(rollAtmFrameCrop(closestAtmPointsInd));
            % The mean aircraft roll associated with the points closest to
            % nadir (decimal degrees).
            
            %..............................................................
            % ATM: Horizontal positions.
            %..............................................................
            distAtmClosest     = dist(closestAtmPointsInd);
            % The vector of distances from nadir associated with the points
            % closest to nadir (metres).
            
            dtsAtmDistPre      = computeSurfaceStats(distAtmClosest,'sipNo');
            % Calculates the descriptive statistics associated with the
            % distance points closest to nadir.
            
            dtsAtmDist{ind}    = dtsAtmDistPre;
            % Put the descriptive statisitcs output structure into the
            % corresponding cell array.
            
            %..............................................................
            % ATM: Vertical positions.
            %..............................................................
            elevationAtmClosest = elevationAtmFrameCrop(closestAtmPointsInd);
            % The vector of surface elevations associated with the points
            % closest to nadir (metres).
            
            dtsAtmElevPre       = computeSurfaceStats(elevationAtmClosest,'sipNo');
            % Calculates the descriptive statistics associated with the
            % surface elevation points closest to nadir.
            
            dtsAtmElev{ind}     = dtsAtmElevPre;
            % Put the descriptive statisitcs output structure into the
            % corresponding cell array.
            
            h_topo(ind) = dtsAtmElevPre.Percentile0950-dtsAtmElevPre.Percentile0050;
            % Calculates the h_topo parameter for each trace in the
            % currently loaded snow radar file
            
            floeElevWav(ind) = nan;
            % Calculates the wavelet derived surface inflection point
            % for each trace in the currently loaded snow radar file
            %..............................................................
            
            %--------------------------------------------------------------
            % If the number of surface elevation points closest to nadir is
            % less than the required number of ATM Points (nominally 200).
            %--------------------------------------------------------------
        elseif length(sortDist)<numAtmPointsForStats
            
            meanAtmPointUtcSecs(ind) = nan;
            meanAtmPointLat(ind)     = nan;
            meanAtmPointLon(ind)     = nan;
            
            meanAtmPitch(ind)        = nan;
            meanAtmRoll(ind)         = nan;
            
            dtsAtmDist{ind}          = nan;
            dtsAtmElev{ind}          = nan;
            
            h_topo(ind)              = nan;
            floeElevWav(ind)         = nan;
            
        end
        
        %------------------------------------------------------------------
        % If there is not positional data for the current trace.
        %------------------------------------------------------------------
    elseif isempty(latAtmFrameCrop)
        
        meanAtmPointUtcSecs(ind) = nan;
        meanAtmPointLat(ind)     = nan;
        meanAtmPointLon(ind)     = nan;
        
        meanAtmPitch(ind)        = nan;
        meanAtmRoll(ind)         = nan;
        
        dtsAtmDist{ind}          = nan;
        dtsAtmElev{ind}          = nan;
        
        h_topo(ind)              = nan;
        floeElevWav(ind)         = nan;
        
    end
    %----------------------------------------------------------------------
    
end

%==========================================================================
%% Outputs.
%==========================================================================

%--------------------------------------------------------------------------
% Input properties.
%--------------------------------------------------------------------------
surfaceAttributes.frameLatitude  = frame.Latitude;
surfaceAttributes.frameLongitude = frame.Longitude;
surfaceAttributes.frameGpsTime   = frame.GPS_time;
surfaceAttributes.frameFastTime  = frame.Time;

surfaceAttributes.radarType      = frame.radarType;
surfaceAttributes.dateStr        = frame.dateStr;
surfaceAttributes.segmentNumber  = frame.segmentNumber;
surfaceAttributes.frameNumber    = frame.frameNumber;
surfaceAttributes.fullFilePath   = frame.fullFilePath;

%--------------------------------------------------------------------------
% Time properties.
%--------------------------------------------------------------------------
surfaceAttributes.frameUtcSecs       = frame.utcSecs;
surfaceAttributes.frameMinUtcSecs    = frame.minUtcSecs;
surfaceAttributes.frameMaxUtcSecs    = frame.maxUtcSecs;
surfaceAttributes.frameMedianUtcSecs = frame.medianUtcSecs;

%--------------------------------------------------------------------------
% Echogram properties.
%--------------------------------------------------------------------------
surfaceAttributes.echogramNumberOfRows = frame.echogramNumberOfRows;
surfaceAttributes.echogramNumberOfCols = frame.echogramNumberOfCols;

%--------------------------------------------------------------------------
% Radar properties.
%--------------------------------------------------------------------------
surfaceAttributes.frameRadarBandwidth = frame.radarBandwidth;

%--------------------------------------------------------------------------
% Surface processing properties.
%--------------------------------------------------------------------------
surfaceAttributes.atmTimeCrop          = atmTimeCrop;
surfaceAttributes.nadirDistanceCrop    = nadirDistanceCrop;
surfaceAttributes.numAtmPointsForStats = numAtmPointsForStats;

%--------------------------------------------------------------------------
% Average surface properties.
%--------------------------------------------------------------------------
surfaceAttributes.meanAtmPointUtcSecs = meanAtmPointUtcSecs;
surfaceAttributes.meanAtmPointLat     = meanAtmPointLat;
surfaceAttributes.meanAtmPointLon     = meanAtmPointLon;
surfaceAttributes.meanAtmPitch        = meanAtmPitch;
surfaceAttributes.meanAtmRoll         = meanAtmRoll;

%--------------------------------------------------------------------------
% Surface descriptive statisitics cell properties.
%--------------------------------------------------------------------------
surfaceAttributes.dtsAtmDist          = dtsAtmDist;
surfaceAttributes.dtsAtmElev          = dtsAtmElev;

%--------------------------------------------------------------------------
% Key surface properties.
%--------------------------------------------------------------------------
surfaceAttributes.h_topo              = h_topo;
surfaceAttributes.floeElevWav         = floeElevWav;

%--------------------------------------------------------------------------
% Processing tracking properties.
%--------------------------------------------------------------------------
surfaceAttributes.dateFrameLoaded     = frame.dateLastLoaded;
surfaceAttributes.dateProcessed       = datestr(clock);
surfaceAttributes.frameTrackingNumber = frame.trackingNumber;

%--------------------------------------------------------------------------
%==========================================================================

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% END OF CODE
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end