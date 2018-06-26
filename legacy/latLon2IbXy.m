function [xMetres,yMetres] = latLon2IbXy(latitude,longitude)
% latLon2IbXy - This function converts (latitude, longitude) into 
% polar stereographic (x,y) coordinates.
%
%--------------------------------------------------------------------------
% Arctic
%--------------------------------------------------------------------------
% The standard projection parameters for Northern hemisphere Operation 
% IceBridge data are:
% 
% Polar Stereographic
% Standard Parallel 70N
% Longitude of the origin (central meridian): 45W
% WGS 84 ellipsoid
% 
% This projection is defined as EPSG:3413.
%--------------------------------------------------------------------------
% Antarctic
%--------------------------------------------------------------------------
% The standard projection parameters for Southern hemisphere Operation 
% IceBridge data are:
% 
% Polar Stereographic
% Standard Parallel 71S
% Longitude of the origin (central meridian): 0
% WGS 84 ellipsoid
% 
% This projection is defined as EPSG:3031. 
% -------------------------------------------------------------------------
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    latitude  - The latitude of the point.
%    longitude - The longitude of the point.
%
% Outputs:
%    xMetres - The polar sterographic x coordinate.
%    yMetres - The polar sterographic y coordinate.
%
% Example: 
% -
%
% Other m-files required: none
%
% Subfunctions: none
%
% MAT-files required:
% - /map/map/almanac.m
% - /map/mapdisp/defaultm.m
% - /map/mapproj/projfwd.m
% - /map/mapproj/stereo.m
% - /nanmean.m
%
% See also: functions/scripts using iceBridge data
%
%
% Author: Thomas Newman, Ph.D., Earth Sciences
%
% February 2015; Last revision: Base. 

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% BEGIN CODE
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%==========================================================================
%% Set up map projection depending if data is in the northern or southern hemisphere
%==========================================================================
if nanmean(latitude) > 0
    % northern hemisphere
    
    mstruct_pole = defaultm('stereo');
    mstruct_pole.origin = [90 -45 0];
    mstruct_pole.geoid = almanac('earth','wgs84','kilometers');
    mstruct_pole.falsenorthing = 0;
    mstruct_pole.falseeasting = 0;
    mstruct_pole.scalefactor = 0.969858730377;
    
    mstruct_pole = defaultm(stereo(mstruct_pole));
    % fix errors if any
    
elseif nanmean(latitude) < 0
    % southern hemisphere
    
    mstruct_pole = defaultm('stereo');
    mstruct_pole.origin = [-90 0 0];
    mstruct_pole.geoid = almanac('earth','wgs84','kilometers');
    mstruct_pole.falsenorthing = 0;
    mstruct_pole.falseeasting = 0;
    mstruct_pole.scalefactor = 0.9727715381;
    
    mstruct_pole = defaultm(stereo(mstruct_pole));
    % fix errors if any
    
end

%==========================================================================
%% calculate projected coordinates
%==========================================================================
[xKm,yKm] = projfwd(mstruct_pole,latitude,longitude);
%converts (latitude, longitude) into polar stereographic (x,y) coordinates.

xMetres = xKm*1000;
yMetres = yKm*1000;
% Convert from kilometers to metres.
%==========================================================================

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% END OF CODE
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end