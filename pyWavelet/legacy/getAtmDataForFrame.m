function atm = getAtmDataForFrame(frame)
% getAtmDataForFrame - get the ATM data corresponding to the medianFrameUtcSecs
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    frame - This is a structure containing the currently loaded frame.
% Outputs:
%    atm - A structure containing the atm points and ancillary information.
%
% Example: 
% -
%
% Other m-files required: 
% - /shared/maputils/wrapTo180.m
% - calcLeapSecond
% - readIlatm1b
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
%% Setup.
%==========================================================================
dateStr = frame.dateStr;
% The date string.

gpsAheadUtc = calcLeapSecond(dateStr);
% The number of seconds the gps time is ahead of the utc time.

numAdjacentFiles = 5;
% The number of adjacent ATM files to concatenate together.

%==========================================================================
%% File indexing.
%==========================================================================
foldersStruct = dir(frame.atmPath);
% List the folder contents as a structure.

foldersCell   = struct2cell(foldersStruct);
% Convert folder structure into a cell array.

nonFileNameCount = 0;
outInd           = 0;
fileNames        = [];
timeIndex        = [];
% Indexing indices. 

%==========================================================================
%% Collecting the nearest ATM files.
%==========================================================================
for ind = 1:size(foldersCell,2)
    
    cellFileNameStr = foldersCell{1,ind};
    % The current filename.
    
    if isempty(strfind(cellFileNameStr,'ILATM1B'))
        
        nonFileNameCount= nonFileNameCount+1;
        % If no .qi or .h5 filename is detected.
        
    elseif ~isempty(strfind(cellFileNameStr,'ILATM1B')) 
        
        outInd              = outInd+1;
        % Advance outInd by one;
        
        fileNames{outInd,1} = cellFileNameStr;
        % Output a cell array of filenames.
        
        startTimeInd = strfind(cellFileNameStr, sprintf('%s_',dateStr))+length(sprintf('%s_',dateStr));
        % The start time index of the filename.
        
        endTimeInd   = strfind(cellFileNameStr, '.ATM')-1;
        % The end time index of the filename.
        
        tempHH = str2num(cellFileNameStr(startTimeInd:startTimeInd+1));
        % The file start hours (from filename).
        
        tempMM = str2num(cellFileNameStr(startTimeInd+2:startTimeInd+3));
        % The file start minutes (from filename).
        
        tempSS = str2num(cellFileNameStr(startTimeInd+4:startTimeInd+5));
        % The file start seconds (from filename).
        
        utcTimeSecs = (tempHH*60*60)+(tempMM*60)+(tempSS) - gpsAheadUtc;
        % Convert GPS time into UTC time.
        
        timeIndex(outInd,1) = utcTimeSecs;
        % The filename DOES NOT correspond to exactly the same start time
        % in the file.
        
    end
end

deltaAtmRadar    = abs(timeIndex - frame.medianUtcSecs);
% Finds the difference between the median frame utc time and the time index
% which represents the approx start time of each of the ATM files.

[deltaAtmRadarVal,deltaAtmRadarInd] = sort(deltaAtmRadar);
% Sorts the differences between the median frame utc time and the time index
% for each of the ATM files.

fileNumbers = sort(deltaAtmRadarInd(1:numAdjacentFiles));
% The filenumber curresponding the to closest (numAdjacentFiles) in time to
% the medianFrameUtcSecs.

%==========================================================================
%% Output the relevant ATM data.
%==========================================================================
latitude      = []; 
longitude     = []; 
elevation     = [];
scanAzimuth   = [];
pitch         = [];
roll          = [];
gpsTimeHHMMSS = [];
% Preallocation.

for ind = 1:length(fileNumbers);

    cellFileNameStr = foldersCell{1, nonFileNameCount+fileNumbers(ind)};
    % the current file.
    
    filePath = strcat(frame.atmPath,cellFileNameStr);
    % The filepath of current file.
    
    atmPre = readIlatm1b(filePath);
    % Read the ATM file.
    
    latitude      = [latitude      , atmPre.latitude];
    longitude     = [longitude     , atmPre.longitude];
    elevation     = [elevation     , atmPre.elevation];
    scanAzimuth   = [scanAzimuth   , atmPre.scanAzimuth];
    pitch         = [pitch         , atmPre.pitch];
    roll          = [roll          , atmPre.roll];
    gpsTimeHHMMSS = [gpsTimeHHMMSS , atmPre.gpsTimeHHMMSS]; 
    % Concatinates the ATM information.
end

longitude  = wrapTo180(longitude);
% Wrap angle in degrees to [-180 180]

%--------------------------------------------------------------------------
%% Converts from gpsTimeHHMMSS to UTC time
%--------------------------------------------------------------------------
MSrem = rem(gpsTimeHHMMSS,1e3);
MS = MSrem;
% Milliseconds. 

SSrem = rem(gpsTimeHHMMSS,1e5);
SS = (SSrem - MSrem)/1e3;
% Seconds. 

MMrem= rem(gpsTimeHHMMSS,1e7);
MM = (MMrem - SSrem)/1e5;
% Minutes. 

HHrem = rem(gpsTimeHHMMSS,1e9);
HH = (HHrem - MMrem)/1e7;
% Hours.

utcTime = ((HH*3600) + (MM*60) + SS +(MS/1000) - gpsAheadUtc)';
% Convert recombine to produce utcTime.
%--------------------------------------------------------------------------

%==========================================================================
%% Construct struct file for output.
%==========================================================================
atm.latitude    = latitude;
atm.longitude   = longitude;
atm.elevation   = elevation;
atm.scanAzimuth = scanAzimuth;
atm.pitch       = pitch;
atm.roll        = roll;
atm.utcTime     = utcTime;
% Output data structure.

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% END OF CODE
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end