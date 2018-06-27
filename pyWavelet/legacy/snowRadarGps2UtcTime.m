function utcTime = snowRadarGps2UtcTime(dateStr,gpsTime)
% snowRadarGps2UtcTime - This function converts the GPS time (seconds) 
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    dateStr - The standard date string [YYYYMMDD]
%    gpsTime - The gps time in seconds.
%
% Outputs:
%    utcTime - The utc time in seconds.
%
% Example: 
% - calcLeapSecond
%
% Other m-files required: 
% 
%
% Subfunctions: none
%
% MAT-files required: none
%
% See also: loadSnowRadarFrame
% 
%
% Author: Thomas Newman, Ph.D., Earth Sciences
%
% February 2015; Last revision: Base. 

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% BEGIN CODE
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%==========================================================================
%% The reference date vect.
%==========================================================================
refDateVect = [1970,1,1,0,0,0];
% The c standard reference time.

%==========================================================================
%% Converting the datestr into numbers.
%==========================================================================
yearNum  = str2num(dateStr(1:4));
monthNum = str2num(dateStr(5:6));
dayNum   = str2num(dateStr(7:8)); 

dateVect = [yearNum,monthNum,dayNum,0,0,0];
% The date vect as specified by the date string.

%==========================================================================
%% Calculating the utc time.
%==========================================================================
gpsAheadUtc = calcLeapSecond(dateStr);
% The number of leap seconds the GPS time is ahead of the UTC time.

elapsedTime = etime(dateVect,refDateVect);
% The time in seconds between the beginning of the day and the reference
% time.

utcTime = gpsTime - elapsedTime - gpsAheadUtc;
% The output utc time in seconds.
%==========================================================================

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% END OF CODE
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end