function gpsAheadUtc = calcLeapSecond(dateStr)
% calcLeapSecond - Calculates the leap second applicable for the current date
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    dateStr - The current date in a date string 
%
% Outputs:
%    gpsAheadUtc - The number of seconds the gps time is ahead of the utc time
%
% Example: 
% -
%
% Other m-files required: none
%
% Subfunctions: none
%
% MAT-files required: none
%
% See also: getAtmDataForFrame snowRadarGps2UtcTime
%
%
% Author: Thomas Newman, Ph.D., Earth Sciences
% 
% February 2015; Last revision: Base.

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% BEGIN CODE
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%==========================================================================
%% Converting the datestr into numbers.
%==========================================================================
yearNum  = str2num(dateStr(1:4));
monthNum = str2num(dateStr(5:6));
dayNum   = str2num(dateStr(7:8));

%==========================================================================
%% Calculating the number of leap seconds to add.
%==========================================================================
if (yearNum >= 2006) && (yearNum < 2009)
    gpsAheadUtc = 14;
    
elseif (yearNum >= 2009) && (yearNum < 2012)
    gpsAheadUtc = 15;
    
elseif (yearNum >= 2012) && (yearNum < 2015)
    
    if yearNum == 2012
        if monthNum >= 7
            gpsAheadUtc = 16;
        elseif monthNum < 7
            gpsAheadUtc = 15;
        end
    elseif yearNum > 2012
        gpsAheadUtc = 16;
    end
    
elseif (yearNum >= 2015)
    
    if yearNum == 2015
        if monthNum >= 7
            gpsAheadUtc = 17;
        elseif monthNum < 7
            gpsAheadUtc = 16;
        end
    elseif yearNum > 2015
        gpsAheadUtc = 17;
    end
    
end

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% END OF CODE
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end