function frame = loadSnowRadarFrame(iceBridgePath, radarType,dateStr,segmentNum,frameNum,loadMode)
% loadSnowRadarFrame - This loads the specified snow radar frame and
% outputs a structure containing the echogram and associated data.
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    radarType  - The snow radar echogram type: 
%                 ('Snow_Radar' or 'Snow_Radar_Deconvolved').
%    dateStr    - The date string of the flight.
%    segmentNum - The segment number of the frame.
%    frameNum   - The frame number.
%
% Outputs:
%    frame - A structure containing the data contained within the frame and
%            ancillary information.
%
% Example: 
% -   
%
% Other m-files required:
% - snowRadarGps2UtcTime
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
% Februrary 2015; Last revision: Base.

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% BEGIN CODE
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%==========================================================================
%% Inputs.
%==========================================================================
segmentStr    = sprintf('%02i', segmentNum);
% Converts the segmentNum into a string using the correct formatting.  

frameStr      = sprintf('%03i', frameNum);
% Converts the frameStr into a string using the correct formatting. 

radarFilePath = sprintf('%s/%s/%s/%s/',iceBridgePath, dateStr(1:4), dateStr, radarType);
% Produces the radar file path pointing to the location of folder.

radarDataPath  = sprintf('%s_%s/Data_%s_%s_%s.mat', dateStr, segmentStr, dateStr, segmentStr, frameStr);
% Produces the file path pointing to the location of data.

fullFilePath   = strcat(radarFilePath, radarDataPath);
% The file path that points to the location of the data.

%==========================================================================
%% loads the specified frame.
%==========================================================================
switch loadMode
    case 'lte'
    frame = load(fullFilePath,'-regexp', '^(?!Data)...');
     % Loads the specified frame.    
        
    case 'full'
     frame = load(fullFilePath);
     % Loads the specified frame.
     
     case 'justData'
     frame = load(fullFilePath,'Data');
     return
end
%==========================================================================
%% Calculate ancillary information.
%==========================================================================
frameUtcSecs         = snowRadarGps2UtcTime(dateStr,frame.GPS_time);
% Converts the frame's GPS times into utc times.

frameMinUtcSecs      = min(frameUtcSecs);
% The minimum utc time in the current frame.

frameMaxUtcSecs      = max(frameUtcSecs);
% The maximum utc time in the current frame.

frameMedianUtcSecs   = median(frameUtcSecs);
% The median utc time in the current frame.

switch loadMode
    case 'full'
        
        echogramNumberOfRows = size(frame.Data,1);
        % The number of rows in echogram.
        
        echogramNumberOfCols = size(frame.Data,2);
        % The number of columns in echogram.
end
%--------------------------------------------------------------------------
% The radar bandwidth.
%--------------------------------------------------------------------------
% switch radarType
%     case 'Snow_Radar'
%         switch dateStr(1:4)
%             case '2011'
%                 switch dateStr
%                     case '20110323'
%                         bandwidth = abs((frame.param_qlook.radar.f1 - ...
%                             frame.param_qlook.radar.f0) * ...
%                             frame.param_qlook.radar.fmult);
%                     otherwise
%                         bandwidth = abs((frame.param_radar.f1 - ...
%                             frame.param_radar.f0) * ...
%                             frame.param_radar.fmult);
%                 end
%         end
%     case 'Snow_Radar_Deconvolved'
%         bandwidth = abs((frame.param_records.radar.wfs.f1-frame.param_records.radar.wfs.f0)*frame.param_records.radar.wfs.fmult);
% end

bandwidth = abs((frame.param_records.radar.wfs.f1-frame.param_records.radar.wfs.f0)*frame.param_records.radar.wfs.fmult);

% The exact bandwidth of the radar system.
%--------------------------------------------------------------------------

%==========================================================================
%% Output.
%==========================================================================

%--------------------------------------------------------------------------
%% Input properties.
%--------------------------------------------------------------------------
frame.radarType     = radarType;
frame.dateStr       = dateStr;
frame.segmentNumber = segmentNum;
frame.frameNumber   = frameNum;
frame.fullFilePath  = fullFilePath;

%--------------------------------------------------------------------------
%% Time properties.
%--------------------------------------------------------------------------
frame.utcSecs       = frameUtcSecs;
frame.minUtcSecs    = frameMinUtcSecs;
frame.maxUtcSecs    = frameMaxUtcSecs;
frame.medianUtcSecs = frameMedianUtcSecs;

%--------------------------------------------------------------------------
%% Radar properties.
%--------------------------------------------------------------------------
frame.radarBandwidth = bandwidth;

%--------------------------------------------------------------------------
%% Echogram properties.
%--------------------------------------------------------------------------
switch loadMode
    case 'full'
        frame.echogramNumberOfRows = echogramNumberOfRows;
        frame.echogramNumberOfCols = echogramNumberOfCols;
end
%--------------------------------------------------------------------------
%% Processing tracking properties.
%--------------------------------------------------------------------------
switch loadMode
    case 'full'
        frame.dateLastLoaded = datestr(clock);
        frame.trackingNumber = round(rand*1e7);
end
%--------------------------------------------------------------------------
%==========================================================================

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% END OF CODE
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end