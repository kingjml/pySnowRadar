% This script loads snow radar frames and picks snow layers using
% wavelet techniques.
% For each picked snow depth associated topographic information is
% calculated that can be used to filter the snow depths.
%
% Example:
% -
%
% Other m-files required:
% - /map/map/almanac.m
% - /map/map/ingeoquad.m
% - /map/map/wgs84Ellipsoid.m
% - /map/mapgeodesy/geodetic2ned.m
% - /map/mapgeodesy/ned2geodetic.m
% - calcSurfaceAttributes
% - generateQcFigureQf
% - loadSnowRadarFrame
% - pickSnowLayer
% - qcUncertCalc
%
% Subfunctions: none
%
% MAT-files required:
% - `resume' files of the form: `resume20110316.mat'
%
% See also: 
%
%
% Author: Thomas Newman, Ph.D., Earth Sciences
%
% February 2015; Last revision: Base.

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% BEGIN CODE
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%==========================================================================
%% Clean up.
%==========================================================================
clc;
myStops = dbstatus('-completenames');
save('myBreakpoints.mat','myStops');
close all;
clear -global;
clear all force;
clear classes;
load('myBreakpoints.mat');
dbstop(myStops);
clear myStops;
delete('myBreakpoints.mat');
% disp('Cleaned up');

addpath(genpath('/Users/tom/Documents/MATLAB'))

%==========================================================================
%% Initialization.
%==========================================================================
colordef white
scrsz      = get(0,'ScreenSize');
wgs84Param = almanac('earth','wgs84','meters');

%==========================================================================
%% File processing limits.
%==========================================================================

%--------------------------------------------------------------------------
% UTC time.
%--------------------------------------------------------------------------
startUtcSecs = [];
% Start UTC time.

endUtcSecs   = [];
% End UTC time.

%--------------------------------------------------------------------------
% Geographic.
%--------------------------------------------------------------------------
minLatPro = []; 
% The minimum latitude to process.

maxLatPro = []; 
% The maximum latitude to process.

minLonPro = []; 
% The minimum longitude to process.

maxLonPro = []; 
% The maximum longitude to process.

boundaryRegionM = [];
% The boundary region around the coordinates.

%==========================================================================
%% CreSiS quality flags.
%==========================================================================
% goodData          = '00000000';
%
% landOrIceberg     = '10000000';
% unclassifiedError = '01000000';
% lowSnr            = '00100000';
% noGoodData        = '00010000';
% missingData       = '00001000';
% verticalStripes   = '00000100';
% deconvolution     = '00000010';
% coherentNoise     = '00000001';

%==========================================================================
%% Resume processing from the last previously processed file.
%==========================================================================
processingMode = 'concat'; % ('picker','qc','figures','concat');
% picker   : Uses wavelet techniques to pick the interfaces.
% qc       : Performs quality control on the picked interfaces.
% concat   : Concatinates the processed files together.
% figures  : Produces quality control figures.
% qualFlag : Generates the quality flag for each flight date.

programMode = 'new'; % ('setup','new','resume','gap');
% Setup  : Makes the raw and processed folders.
% New    : Processes from scratch.
% Resume : Resume processing from last processed file.
% Gap    : Processes any files not yet processed.

disp(sprintf('Program mode: %s',programMode))

switch programMode
    case 'setup'
        str = input('This will make icebridge\n folder structure, continue (y/n): ','s');
    case 'new'
        str = input('This will overwrite the previously\nsaved files, continue (y/n): ','s');
    case 'resume'
        str = input('This will resume processing\nfrom the last processed file, continue (y/n): ','s');
    case 'gap'
        str = input('This fill in any gaps in processing\nonly process absent files, continue (y/n): ','s');
end

if strcmp(str,'y') == 0
    disp('Program ending')
    return
end

%==========================================================================
%% Define the file paths.
%==========================================================================
radarType     = 'Snow_Radar'; 
% Snow_Radar_Deconvolved = picks interfaces.
% Snow_Radar             = Does now pick interfaces.
 
iceBridgePath = '/Volumes/tnDataDisk/iceBridgeData';

%--------------------------------------------------------------------------
% 2009.
%--------------------------------------------------------------------------
% dateStrVect2009 = {...
%     '20090331',...
%     '20090402',...
%     '20090405',...
%     '20090421',...
%     '20090425',...
%     };

%--------------------------------------------------------------------------
% 2010.
%--------------------------------------------------------------------------
% dateStrVect2010 = {...
%     '20100323',...
%     '20100326',...
%     '20100402',...
%     '20100405',...
%     '20100412',...
%     '20100419',...
%     '20100420',...
%     '20100421',...
%     };

%--------------------------------------------------------------------------
% 2011.
%--------------------------------------------------------------------------
% dateStrVect2011 = {...
%     '20110316',...
%     '20110317',...
%     '20110318',...
%     '20110322',...
%     '20110323',...
%     '20110325',...
%     '20110326',...
%     '20110328',...
%     '20110415',...
%     };

%--------------------------------------------------------------------------
% 2012.
%--------------------------------------------------------------------------
% dateStrVect2012 = {...
%     '20120314',...
%     '20120315',...
%     '20120316',...
%     '20120317',...
%     '20120319',...
%     '20120321',...
%     '20120322',...
%     '20120323',...
%     '20120326',...
%     '20120327',...
%     '20120328',...
%     '20120329',...
%     '20120402',...
%     };

dateStrVect2012 = {...
    '20120410',...
    };

%--------------------------------------------------------------------------
% 2013.
%--------------------------------------------------------------------------
dateStrVect2013 = {...
    '20130321',...
    '20130322',...
    '20130323',...
    '20130324',...
    '20130326',...
    '20130327',...
    '20130422',...
    '20130424',...
    '20130425',...
    };

%--------------------------------------------------------------------------
% 2014.
%--------------------------------------------------------------------------
% dateStrVect2014 = {...
%     '20140313',...
%     '20140315',...
%     '20140317',...
%     '20140318',...
%     '20140319',...
%     '20140321',...
%     '20140324',...
%     '20140325',...
%     '20140326',...
%     '20140328',...
%     '20140331',...
%     '20140403',...
%     '20140428',...
%     };

dateStrVect2014 = {...
    '20140312',...
    '20140314',...
    };

%--------------------------------------------------------------------------
% 2015.
%--------------------------------------------------------------------------
% dateStrVect2015 = {...
%     '20150319',...
%     '20150324',...
%     '20150325',...
%     '20150326',...
%     '20150327',...
%     '20150329',...
%     '20150330',...
%     '20150401',...
%     '20150403',...
%     '20150508',...
%     };
%--------------------------------------------------------------------------

dateStrVect =  [...
    dateStrVect2013,...
    dateStrVect2014,...
    dateStrVect2012];


%dateStrVect = dateStrVect2014;
%==========================================================================
%% Block process the snow radar data.
%==========================================================================
for flightInd = 1:length(dateStrVect);
    % For each flight.
    
    dateStr = dateStrVect{flightInd};
    % The date string
    
    %----------------------------------------------------------------------
    % Load the resume .mat file.
    %----------------------------------------------------------------------
    clear resume
    
    resumeFolderPath = sprintf('%s/%s/%s/',...
        iceBridgePath, dateStr(1:4),'resumeFiles');
    
    switch processingMode
        case 'picker'
            resumeFilePath = sprintf('%sresume%s.mat',...
                resumeFolderPath,dateStr);
            
        case 'qc'
            resumeFilePath = sprintf('%sresumeQc%s.mat',...
                resumeFolderPath,dateStr);
            
        case 'figures'
            resumeFilePath = sprintf('%sresumeFig%s.mat',...
                resumeFolderPath,dateStr);
            
        case 'concat'
            resumeFilePath = sprintf('%sresumeConcat%s.mat',...
                resumeFolderPath,dateStr);
            
        case 'qualFlag'
            resumeFilePath = sprintf('%sresumeQualFlag%s.mat',...
                resumeFolderPath,dateStr);
            
    end
    
    switch programMode
        case 'setup'
            if exist(resumeFolderPath) == 0
                mkdir(resumeFolderPath)
            end
            
        case {'new','gap'}
            segmentIndSt = 1;
            frameIndSt   = 1;
            
        case 'resume'
            if exist(resumeFilePath) == 2
                load(resumeFilePath)
                segmentIndSt = resume.segmentInd;
                frameIndSt   = resume.frameInd;
                
            elseif exist(resumeFilePath) ~= 2
                segmentIndSt = 1;
                frameIndSt   = 1;
            end
    end
    
    %----------------------------------------------------------------------
    %% Segment directory.
    %----------------------------------------------------------------------
    segmentPathPre = sprintf('%s/%s/%s/',...
        iceBridgePath,...
        dateStr(1:4), ...
        dateStr);
    % The file path to the segment folders.
    
    switch programMode
        case 'setup'
            if (exist(segmentPathPre)==0)
                
                segmentPath = sprintf('%s%s',segmentPathPre,radarType);
                mkdir(segmentPath)
                
                atmPath     = sprintf('%s%s',segmentPathPre,'ATM');
                mkdir(atmPath)
                
                qualityFlagPath = sprintf('%s%s',segmentPathPre,'Quality_Flag');
                mkdir(qualityFlagPath)
                
            end
            
        case {'new','gap','resume'}
            segmentPath = sprintf('%s%s',segmentPathPre,radarType);
            segmentsDir = dir(segmentPath);
            % Lists the files and folders in the segmentPath directory.
            
            atmPath         = sprintf('%s%s/',segmentPathPre,'ATM');
            qualityFlagPath = sprintf('%s%s/',segmentPathPre,'Quality_Flag');
    end
    
    switch programMode
        case 'setup'
            continue
    end
    
    %----------------------------------------------------------------------
    % Preallocation for concat data.
    %----------------------------------------------------------------------
    switch processingMode
        case 'concat'
            
            %..............................................................
            % Book keeping.
            %..............................................................
            dateNumVect                 = [];
            segmentNumVect              = [];
            frameNumVect                = [];
            traceNumVect                = [];
            utcTimeVect                 = [];
            latitudeVect                = [];
            longitudeVect               = [];
            
            %..............................................................
            % Quality control.
            %..............................................................
            qcArray                     = [];
            atmDataAvailableVect        = [];
            snowDataAvailableVect       = [];
            
            % QC2..........................................................
            minDistAtmPoint2NadirVect   = [];
            minDistAtmPoint2NadirQcVect = [];
            
            % QC3..........................................................
            maxDistAtmPoint2NadirVect   = [];
            maxDistAtmPoint2NadirQcVect = [];
            
            % QC4..........................................................
            meanAtmPitchVect            = [];
            meanAtmPitchQcVect          = [];
            
            % QC5..........................................................
            meanAtmRollVect             = [];
            meanAtmRollQcVect           = [];
            
            % QC6..........................................................
            hTopoVect                   = [];
            hTopoQcVect                 = [];
            
            % QC7..........................................................
            interfaceOrderVect          = [];
            interfaceOrderQcVect        = [];
            
            % QC8..........................................................
            minSnowDepthQcVect          = [];
            
            % QC9..........................................................
            maxSnowDepthQcVect          = [];
            
            % QC10.........................................................
            snrAirSnowVect              = [];
            minSnrAirSnowQcVect         = [];
            
            % QC11.........................................................
            snrSnowIceVect              = [];
            minSnrSnowIceQcVect         = [];
            
            % QC12.........................................................
            noiseAndClutterPowerVect    = [];
            noisePowerThresholdVect     = [];
            
            %..............................................................
            % Snow depth.
            %..............................................................
            rangeBinAirSnowVect         = [];
            rangeBinSnowIceVect         = [];
            snowDepthVect               = [];
            radarPrecisionVect          = [];
            
            %..............................................................
            % Surface elevation.
            %..............................................................
            elevationMinVect            = [];
            elevation05Vect             = [];
            elevation50Vect             = [];
            elevation95Vect             = [];
            elevationMaxVect            = [];
            elevationMeanVect           = [];
            elevationStdVect            = [];
            numAtmPointsVect            = [];
            
            %..............................................................
            % For intercomparison.
            %..............................................................
            rangeBetweenInterfacesVect  = [];
            rangeToAirSnowInterfaceVect = [];
            rangeToSnowIceInterfaceVect = [];
            fastTimeRangeVect           = [];
            
            %..............................................................
            % Cresis quality flag.
            %..............................................................
            qualityFlagArrayStr = {'landOrIceberg', 'unclassifiedError',...
                'lowSnr', 'noGoodData','missingData','verticalStripes',...
                'deconvolution','coherentNoise'};
            % The column headers for the quality flag array.
            %..............................................................
            
    end
    %----------------------------------------------------------------------
    % Preallocation for quality flag array.
    %----------------------------------------------------------------------
    qualityFlagArray = [];
    % The quality flag array.
    
    qualityFlagInd = 0;
    % The quality flag index.
    
    %----------------------------------------------------------------------
    %% For each segment.
    %----------------------------------------------------------------------
    for segmentInd = segmentIndSt:length(segmentsDir)
        
        if ~isempty(strfind(segmentsDir(segmentInd).name,dateStr))
            % If the segment name matches the segment file format.
            
            %..............................................................
            % Raw data path.
            %..............................................................
            % /Volumes/polarDrive/iceBridgeData/2012/20120314/Snow_Radar_Deconvolved/20120314_03
            
            framePath = sprintf('%s/%s/%s/%s/%s/%s/',...
                iceBridgePath,dateStr(1:4),...
                dateStr,...
                radarType,...
                segmentsDir(segmentInd).name);
            % The file path to the current frames.
            
            framesDir = dir(framePath);
            % Lists the files in the framePath directory.
            % {'Data_20120315_03_001.mat'}
            
            %..............................................................
            % Processed data path.
            %..............................................................
            % /Volumes/polarDrive/iceBridgeData/2012Pro/20120314/Snow_Radar_Deconvolved/20120314_02
            
            proFramePath = sprintf('%s/%sPro/%s/%s/%s/%s/',...
                iceBridgePath,dateStr(1:4),...
                dateStr,...
                radarType,...
                segmentsDir(segmentInd).name);
            % The file path to the current frames.
            
            if (exist(proFramePath)==0)
                mkdir(proFramePath)
                % Makes the processed folder.
                
            end
            
            %..............................................................
            % The figure path.
            %..............................................................
            % /Volumes/polarDrive/iceBridgeData/2012Pro/20120314/Snow_Radar_Deconvolved/20120314_02
            
            figFramePath = sprintf('%s/%sFig/%s/%s/%s/%s/',...
                iceBridgePath,dateStr(1:4),...
                dateStr,...
                radarType,...
                segmentsDir(segmentInd).name);
            % The file path to the current frames.
            
            if (exist(figFramePath)==0)
                mkdir(figFramePath)
                % Makes the processed folder.
                
            end
            
            %..............................................................
            % Load quality control file.
            %..............................................................
            qualityFlagFile = sprintf('frames_%s.mat',segmentsDir(segmentInd).name);
            % The quality control file.
            
            switch radarType
                case 'Snow_Radar_Deconvolved'
                    load(sprintf('%s%s',qualityFlagPath,qualityFlagFile))
                    % Load the quality control data.
            end
            %..............................................................
            
            
            %--------------------------------------------------------------
            % for each frame.
            %--------------------------------------------------------------
            for frameInd = frameIndSt:length(framesDir)
                
                if ~isempty(strfind(framesDir(frameInd).name,'.mat'))
                    % If the frame name matches the frame file format.
                    
                    %......................................................
                    % Raw file name.
                    %......................................................
                    fileName   = framesDir(frameInd).name;
                    % The filename of the frame.
                    % Data_20120315_03_033.mat
                    
                    %......................................................
                    % Processed file name.
                    %......................................................
                    proFileName   = sprintf('proDecon_%s',fileName);
                    % The filename of the frame.
                    
                    wholeProFilePath = strcat(proFramePath,proFileName);
                    % The whole file name.
                    
                    %......................................................
                    % Segment and frame number.
                    %......................................................
                    segmentNum = str2num(fileName(15:16));
                    % The current segment number.
                    
                    frameNum   = str2num(fileName(18:20));
                    % The current frame number.
                    
                    disp(sprintf('date:%s | segment:%2.0f | frame:%3.0f',...
                        dateStr,segmentNum,frameNum));
                    % Displays the current file being processed.
                    
                    %......................................................
                    % Quality flag.
                    %......................................................
                    switch radarType
                        case 'Snow_Radar_Deconvolved'
                            
                            frameQualityFlagDec = frames.quality(frameNum);
                            % The frame quality decimal.
                            
                            frameQualityFlagBinStr = dec2bin(frameQualityFlagDec,8);
                            % The frame quality metric
                            
                        case 'Snow_Radar'
                            frameQualityFlagBinStr = '00000010';
                            % Assign deconvolution flag to non-deconvolved
                            % data.
                    end
                    
                    % goodData          = '00000000';
                    %
                    % landOrIceberg     = '10000000';
                    % unclassifiedError = '01000000';
                    % lowSnr            = '00100000';
                    % noGoodData        = '00010000';
                    % missingData       = '00001000';
                    % verticalStripes   = '00000100';
                    % deconvolution     = '00000010';
                    % coherentNoise     = '00000001';
                    
                    % .....................................................
                    % Flagged as good when actually bad.
                    % .....................................................
                    % 20100326: seg1, file765  , missingData
                    % 20100402: seg1, file22   , missingData
                    % 20100421: seg3, file 140 , missingData
                    % .....................................................
              
                    switch processingMode
                        case 'qualFlag'
                            qualityFlagInd =  qualityFlagInd+1;
                            % The quality flag index.
                            
                            qualityFlagArray(qualityFlagInd,1)  = segmentNum;
                            qualityFlagArray(qualityFlagInd,2)  = frameNum;
                            
                            qualityFlagArray(qualityFlagInd,3)  = str2num(frameQualityFlagBinStr(1));
                            qualityFlagArray(qualityFlagInd,4)  = str2num(frameQualityFlagBinStr(2));
                            qualityFlagArray(qualityFlagInd,5)  = str2num(frameQualityFlagBinStr(3));
                            qualityFlagArray(qualityFlagInd,6)  = str2num(frameQualityFlagBinStr(4));
                            qualityFlagArray(qualityFlagInd,7)  = str2num(frameQualityFlagBinStr(5));
                            qualityFlagArray(qualityFlagInd,8)  = str2num(frameQualityFlagBinStr(6));
                            qualityFlagArray(qualityFlagInd,9)  = str2num(frameQualityFlagBinStr(7));
                            qualityFlagArray(qualityFlagInd,10) = str2num(frameQualityFlagBinStr(8));
                            % The different columns.
                            
                            %..............................................
                            % Amend qc flag.
                            %..............................................
                            switch dateStr
                                case 20100326
                                    if segmentNum == 1
                                        if frameNum == 765
                                            qualityFlagArray(qualityFlagInd,7) = 1;
                                        end
                                    end
                                case 20100402
                                    if segmentNum == 1
                                        if frameNum == 22
                                            qualityFlagArray(qualityFlagInd,7) = 1;
                                        end
                                    end
                                case 20100421
                                    if segmentNum == 3
                                        if frameNum == 140
                                            qualityFlagArray(qualityFlagInd,7) = 1;
                                        end
                                    end
                            end
                            %..............................................
                            
                            
                        case 'concat'
                            
                            qualityFlagVect     = nan(1,8);
                            
                            qualityFlagVect(1)  = str2num(frameQualityFlagBinStr(1));
                            qualityFlagVect(2)  = str2num(frameQualityFlagBinStr(2));
                            qualityFlagVect(3)  = str2num(frameQualityFlagBinStr(3));
                            qualityFlagVect(4)  = str2num(frameQualityFlagBinStr(4));
                            qualityFlagVect(5)  = str2num(frameQualityFlagBinStr(5));
                            qualityFlagVect(6)  = str2num(frameQualityFlagBinStr(6));
                            qualityFlagVect(7)  = str2num(frameQualityFlagBinStr(7));
                            qualityFlagVect(8)  = str2num(frameQualityFlagBinStr(8));
                            % The different columns.
                            
                            %..............................................
                            % Amend qc flag.
                            %..............................................
                            switch dateStr
                                case '20100326'
                                    if segmentNum == 1
                                        if frameNum == 765
                                            qualityFlagVect(5) = 1;
                                            disp('Ammending qc flag')
                                        end
                                    end
                                case '20100402'
                                    if segmentNum == 1
                                        if frameNum == 22
                                            qualityFlagVect(5) = 1;
                                            disp('Ammending qc flag')
                                        end
                                    end
                                case '20100421'
                                    if segmentNum == 3
                                        if frameNum == 140
                                            qualityFlagVect(5) = 1;
                                            disp('Ammending qc flag')
                                        end
                                    end
                            end
                            %..............................................
                            
                    end
                    
                    %......................................................
                    % Save the resume .mat file.
                    %......................................................
                    resume.dateStr    = dateStr;
                    resume.segmentNum = segmentNum;
                    resume.frameNum   = frameNum;
                    
                    resume.segmentInd = segmentInd;
                    resume.frameInd   = frameInd;
                    
                    save(resumeFilePath,'resume')
                    %......................................................
                    
                    
                    %------------------------------------------------------
                    % Main snow radar processing scripts.
                    %------------------------------------------------------
                    switch processingMode
                        
                        %..................................................
                        case {'picker','qc'}
                            %..............................................
                            pathFileType = sprintf('%s/%sPro/%s/%s/%s/%s/',...
                                iceBridgePath,dateStr(1:4),...
                                dateStr,...
                                radarType,...
                                segmentsDir(segmentInd).name);
                            
                            fileNameFileType   = sprintf('proDecon_%s',fileName);
                            % The filename of the frame.
                            
                            wholeFilePathFileType = strcat(pathFileType,fileNameFileType);
                            % The whole file name.
                            
                            %..............................................
                        case 'concat'
                            %..............................................
                            pathFileType = sprintf('%s/%sPro/',...
                                iceBridgePath,dateStr(1:4));
                            
                            fileNameFileType   = sprintf('qualityControlData%s',dateStr);
                            % The filename of the frame.
                            
                            wholeFilePathFileType = strcat(pathFileType,fileNameFileType);
                            % The whole file name.
                            
                            %..............................................
                        case 'figures'
                            %..............................................
                            pathFileType = sprintf('%s/%sFig/%s/%s/%s/%s/',...
                                iceBridgePath,dateStr(1:4),...
                                dateStr,...
                                radarType,...
                                segmentsDir(segmentInd).name);
                            
                            fileNameFileType   = sprintf('Snow_Radar_Deconvolved_%s.png',fileName(6:20));
                            % The filename of the frame.
                            
                            wholeFilePathFileType = strcat(pathFileType,fileNameFileType);
                            % The whole file name.
                            
                            %..............................................
                        case 'qualFlag'
                            %..............................................
                            pathFileType = sprintf('%s/%sPro/',...
                                iceBridgePath,dateStr(1:4));
                            
                            fileNameFileType   = sprintf('qualityFlagArray%s.mat',dateStr);
                            % The filename of the frame.
                            
                            wholeFilePathFileType = strcat(pathFileType,fileNameFileType);
                            % The whole file name.
                            
                            %..............................................
                    end
                    
                    % picker   : Uses wavelet techniques to pick the interfaces.
                    % qc       : Performs quality control on the picked interfaces.
                    % concat   : Concatinates the processed files together.
                    % figures  : Produces quality control figures.
                    % qualFlag : Generates the quality flag for each flight date.
                    
                    %......................................................
                    % Try to process the current frame and if there are
                    % any errors then go to the next frame.
                    %......................................................
                    switch programMode
                        case 'gap'
                            
                            processFile = exist(wholeFilePathFileType)==0;
                        case {'new','resume'}
                            processFile = true;
                    end
                    
                    if processFile==1
                        try
                            
                            switch processingMode
                                case {'picker','qc'}
                                    
                                    %--------------------------------------------------
                                    % Are quality control metrics met?
                                    %--------------------------------------------------
                                    frameLte = loadSnowRadarFrame(iceBridgePath,radarType,dateStr,segmentNum,frameNum,'lte');
                                    % Loads frame in 'lte' mode.
                                    
                                    %..............................................
                                    %% Test to see if within utc time window.
                                    %..............................................
                                    if (~isempty(startUtcSecs)) && (~isempty(endUtcSecs))
                                        inUtcTimeWindow = (frameLte.minUtcSecs>=startUtcSecs) &...
                                            (frameLte.maxUtcSecs<=endUtcSecs);
                                        % Are the utc frame boundaries within the
                                        % specified limits
                                        
                                    elseif (isempty(startUtcSecs)) && (isempty(endUtcSecs))
                                        
                                        inUtcTimeWindow = true;
                                        
                                    end
                                    
                                    %..............................................
                                    %% Test to see if within geographic limits.
                                    %..............................................
                                    if (~isempty(minLatPro)) && (~isempty(maxLatPro))...
                                            && (~isempty(minLonPro)) && (~isempty(maxLonPro))
                                        
                                        lat0 = median([minLatPro,maxLatPro]);
                                        lon0 = median([minLonPro,maxLonPro]);
                                        h0   = 0;
                                        % local coordiante system origin.
                                        
                                        [xNorth,yEast,zDown] = geodetic2ned([minLatPro,maxLatPro],...
                                            [minLonPro,maxLonPro],...
                                            [0, 0],...
                                            lat0,lon0,h0,...
                                            wgs84Ellipsoid);
                                        % Geodetic to local Cartesian NED
                                        
                                        xNorthBound = [min(xNorth)-boundaryRegionM,max(xNorth)+boundaryRegionM];
                                        yEastBound  = [min(yEast)-boundaryRegionM,max(yEast)+boundaryRegionM];
                                        % Add boundary region.
                                        
                                        [latBound,lonBound,hBound] = ned2geodetic(xNorthBound,yEastBound,zDown,lat0,lon0,h0,wgs84Ellipsoid);
                                        % Local Cartesian NED to geodetic
                                        
                                        tfGeo = ingeoquad(frameLte.Latitude, frameLte.Longitude,...
                                            latBound,lonBound);
                                        % True for points inside or on lat-lon quadrangle
                                        
                                        inGeoWindow = any(tfGeo);
                                        % Is any part of the frame within the specified
                                        % geographic limits.
                                        
                                    elseif (isempty(minLatPro)) && (isempty(maxLatPro))...
                                            && (isempty(minLonPro)) && (isempty(maxLonPro))
                                        
                                        inGeoWindow = true;
                                        
                                    end
                                    
                                    %..............................................
                                    % Test to see if data is 'good'.
                                    %..............................................
                                    % goodData          = '00000000';
                                    %
                                    % landOrIceberg     = '10000000';
                                    % unclassifiedError = '01000000';
                                    % lowSnr            = '00100000';
                                    % noGoodData        = '00010000';
                                    % missingData       = '00001000';
                                    % verticalStripes   = '00000100';
                                    % deconvolution     = '00000010';
                                    % coherentNoise     = '00000001';
                                    
                                    isDataQualityGood = str2double(frameQualityFlagBinStr)==0;
                                    % is the data quality good.
                                    
                                    %..............................................
                                    % Combine the flags
                                    %..............................................
                                    processFileTf = inUtcTimeWindow & inGeoWindow; % & isDataQualityGood;
                                    % Process file if contraints are met
                                    %..............................................
                                    
                                case {'concat','figures'}
                                    
                                    processFileTf = true;
                                    % Process file if contraints are met
                                    
                                case 'qualFlag'
                                    processFileTf = false;
                                    
                            end
                            
                            if processFileTf==1
                                
                                switch processingMode
                                    case 'picker'
                                        %--------------------------------------------------
                                        % Load the frame and ancillary information.
                                        disp('...Loading frame...')
                                        %--------------------------------------------------
                                        frame = loadSnowRadarFrame(iceBridgePath,radarType,dateStr,segmentNum,frameNum,'full');
                                        % 'full','lte','justData'
                                        
                                        frame.qualityFlag = frameQualityFlagBinStr;
                                        % The quality metric.
                                        
                                        frame.atmPath         = atmPath;
                                        % The ATM path
                                        
                                        frame.qualityFlagPath = qualityFlagPath;
                                        % The quakity flag path.
                                        
                                        %--------------------------------------------------
                                        % Calculates the surface attributes using the
                                        % information supplied in the frame structure.
                                        disp('...Calculating surface attributes...')
                                        %--------------------------------------------------
                                        surfaceAttributes = calcSurfaceAttributes(frame);
                                        % Calculate the surface attributes.
                                        
                                        %--------------------------------------------------
                                        % Picks the range bins corresponding to the
                                        % air/snow and snow/ice interfaces for each trace
                                        % in the frame.
                                        disp('...Picking snow layer...')
                                        %--------------------------------------------------
                                        snow = pickSnowLayer(frame);
                                        % Picks the snow layer interfaces.
                                        
                                        %--------------------------------------------------
                                        % Saving the processed data.
                                        disp('...Saving file...')
                                        %--------------------------------------------------
                                        frame.Data   = [];
                                        % Removes the echogram data to save space.
                                        
                                        save(wholeProFilePath,'frame','surfaceAttributes','snow')
                                        % Saving the frame, surfaceAttributes and snow
                                        % structures.
                                        clc
                                        %--------------------------------------------------
                                        
                                    case 'qc'
                                        
                                        %--------------------------------------------------
                                        % loading the processed data.
                                        disp('...Loading file...')
                                        %--------------------------------------------------
                                        load(wholeProFilePath,'frame','surfaceAttributes','snow')
                                        % Saving the frame, surfaceAttributes and snow
                                        % structures.
                                        
                                        frameData   = loadSnowRadarFrame(iceBridgePath,radarType,dateStr,segmentNum,frameNum,'justData');
                                        % Loads just the frame data.
                                        
                                        frame.Data   = frameData.Data;
                                        % Adds the echogram data to the frame struct.
                                        
                                        clear frameData
                                        
                                        %--------------------------------------------------
                                        % Performing snow radar data
                                        % quality control
                                        disp('...Performing quality control...')
                                        %--------------------------------------------------
                                        snowQc = qcUncertCalc(frame,surfaceAttributes,snow);
                                        % Performs snow radar quality control
                                        
                                        frame.Data   = [];
                                        % Removes the echogram data to save space.
                                        
                                        save(wholeProFilePath,'frame','surfaceAttributes','snow','snowQc')
                                        % Saving the snowQc struct.
                                        %--------------------------------------------------
                                        
                                    case 'concat'
                                        
                                        %--------------------------------------------------
                                        % Concatinates the snow qc files
                                        disp('...Loading file...')
                                        %--------------------------------------------------
                                        load(wholeProFilePath,'snowQc')
                                        
                                        onesColumn = ones(snowQc.echogramNumberOfCols,1);
                                        % A column of ones.
                                        
                                        %..............................................................
                                        % Book keeping.
                                        %..............................................................
                                        dateNumVect                 = [dateNumVect                 ; snowQc.dateNumVect(:)];
                                        segmentNumVect              = [segmentNumVect              ; snowQc.segmentNumVect(:)];
                                        frameNumVect                = [frameNumVect                ; snowQc.frameNumVect(:)];
                                        traceNumVect                = [traceNumVect                ; snowQc.traceNumVect(:)];
                                        utcTimeVect                 = [utcTimeVect                 ; snowQc.frameUtcSecs(:)];
                                        latitudeVect                = [latitudeVect                ; snowQc.frameLatitude(:)];
                                        longitudeVect               = [longitudeVect               ; snowQc.frameLongitude(:)];
                                        
                                        %..............................................................
                                        % Quality control.
                                        %..............................................................
                                        qcArray                     = [qcArray                     ; snowQc.qcArray.'];
                                        atmDataAvailableVect        = [atmDataAvailableVect        ; snowQc.atmDataAvailable(:)*1];
                                        snowDataAvailableVect       = [snowDataAvailableVect       ; snowQc.snowDataAvailable(:)*1];
                                        
                                        %%QC2..........................................................
                                        minDistAtmPoint2NadirVect   = [minDistAtmPoint2NadirVect   ; snowQc.minDistAtmPoint2NadirVect(:)];
                                        minDistAtmPoint2NadirQcVect = [minDistAtmPoint2NadirQcVect ; onesColumn*snowQc.minDistAtmPoint2NadirQc];
                                        
                                        % QC3..........................................................
                                        maxDistAtmPoint2NadirVect   = [maxDistAtmPoint2NadirVect   ; snowQc.maxDistAtmPoint2NadirVect(:)];
                                        maxDistAtmPoint2NadirQcVect = [maxDistAtmPoint2NadirQcVect ; onesColumn*snowQc.maxDistAtmPoint2NadirQc];
                                        
                                        % QC4..........................................................
                                        meanAtmPitchVect            = [meanAtmPitchVect            ; snowQc.meanAtmPitchVect(:)];
                                        meanAtmPitchQcVect          = [meanAtmPitchQcVect          ; onesColumn*snowQc.meanAtmPitchQc];
                                        
                                        % QC5..........................................................
                                        meanAtmRollVect             = [meanAtmRollVect             ; snowQc.meanAtmRollVect(:)];
                                        meanAtmRollQcVect           = [meanAtmRollQcVect           ; onesColumn*snowQc.meanAtmRollQc];
                                        
                                        % QC6..........................................................
                                        snowQc.hTopoQc = 0.5;
                                        hTopoVect                   = [hTopoVect                   ; snowQc.h_topoVect(:)];
                                        hTopoQcVect                 = [hTopoQcVect                 ; onesColumn*snowQc.hTopoQc];
                                        
                                        % QC7..........................................................
                                        snowQc.interfaceOrderQc     = 1;
                                        interfaceOrderVect          = [interfaceOrderVect          ; snowQc.interfaceOrderVect(:)];
                                        interfaceOrderQcVect        = [interfaceOrderQcVect        ; onesColumn*snowQc.interfaceOrderQc];
                                        
                                        % QC8..........................................................
                                        minSnowDepthQcVect          = [minSnowDepthQcVect          ; onesColumn*snowQc.minSnowDepthQc];
                                        
                                        % QC9..........................................................
                                        maxSnowDepthQcVect          = [maxSnowDepthQcVect          ; onesColumn*snowQc.maxSnowDepthQc];
                                        
                                        % QC10..........................................................
                                        snrAirSnowVect              = [snrAirSnowVect              ; snowQc.snrAirSnowVect(:)];
                                        minSnrAirSnowQcVect         = [minSnrAirSnowQcVect         ; onesColumn*snowQc.minSnrAirSnowQc];
                                        
                                        % QC11..........................................................
                                        snrSnowIceVect              = [snrSnowIceVect              ; snowQc.snrSnowIceVect(:)];
                                        minSnrSnowIceQcVect         = [minSnrSnowIceQcVect         ; onesColumn*snowQc.minSnrSnowIceQc];
                                        
                                        % QC12..........................................................
                                        noiseAndClutterPowerVect    = [noiseAndClutterPowerVect    ; snowQc.noiseAndClutterPowerVect(:)];
                                        noisePowerThresholdVect     = [noisePowerThresholdVect     ; snowQc.noisePowerThresholdVect(:)];
                                        
                                        %..............................................................
                                        % Snow depth.
                                        %..............................................................
                                        rangeBinAirSnowVect         = [rangeBinAirSnowVect         ; snowQc.rangeBinVectAirSnow(:)];
                                        rangeBinSnowIceVect         = [rangeBinSnowIceVect         ; snowQc.rangeBinVectSnowIce(:)];
                                        snowDepthVect               = [snowDepthVect               ; snowQc.snowDepth(:)];
                                        radarPrecisionVect          = [radarPrecisionVect          ; snowQc.radarPrecision(:)];
                                        
                                        %..............................................................
                                        % Surface elevation.
                                        %..............................................................
                                        elevationMinVect            = [elevationMinVect            ; snowQc.elevationMinVect(:)];
                                        elevation05Vect             = [elevation05Vect             ; snowQc.elevation05Vect(:)];
                                        elevation50Vect             = [elevation50Vect             ; snowQc.elevation50Vect(:)];
                                        elevation95Vect             = [elevation95Vect             ; snowQc.elevation95Vect(:)];
                                        elevationMaxVect            = [elevationMaxVect            ; snowQc.elevationMaxVect(:)];
                                        elevationMeanVect           = [elevationMeanVect           ; snowQc.elevationMeanVect(:)];
                                        elevationStdVect            = [elevationStdVect            ; snowQc.elevationStdVect(:)];
                                        numAtmPointsVect            = [numAtmPointsVect            ; snowQc.numAtmPointsVect(:)];
                                        
                                        %..............................................................
                                        % For intercomparison.
                                        %..............................................................
                                        rangeBetweenInterfacesVect  = [rangeBetweenInterfacesVect  ; snowQc.rangeBetweenInterfaces(:)];
                                        rangeToAirSnowInterfaceVect = [rangeToAirSnowInterfaceVect ; snowQc.rangeToAirSnowInterface(:)];
                                        rangeToSnowIceInterfaceVect = [rangeToSnowIceInterfaceVect ; snowQc.rangeToSnowIceInterface(:)];
                                        fastTimeRangeVect           = [fastTimeRangeVect           ; snowQc.fastTimeRange(:)];
                                        
                                        qualityFlagArray            = [qualityFlagArray            ; repmat(qualityFlagVect,snowQc.echogramNumberOfCols,1)];
                                        %..............................................................
                                        
                                        %--------------------------------------------------
                                        
                                    case 'figures'
                                        
                                        %--------------------------------------------------
                                        % loading the processed data.
                                        disp('...Generating figure...')
                                        %--------------------------------------------------
                                        load(wholeProFilePath,'frame','surfaceAttributes','snow','snowQc')
                                        % Saving the frame, surfaceAttributes and snow
                                        % structures.
                                        
                                        frameData   = loadSnowRadarFrame(iceBridgePath,radarType,dateStr,segmentNum,frameNum,'justData');
                                        % Loads just the frame data.
                                        
                                        frame.Data   = frameData.Data;
                                        % Adds the echogram data to the frame struct.
                                        
                                        clear frameData
                                        
                                        %outputFigurePath = generateQcFigure(frame,surfaceAttributes,snow,snowQc,figFramePath);
                                        outputFigurePath = generateQcFigureQf(frame,surfaceAttributes,snow,snowQc,figFramePath);
                                        % Generates the qc figure.
                                        %--------------------------------------------------
                                        
                                end
                                
                                
                            elseif processFileTf==0
                                disp('Not processed');
                                % Save error message as a txt file
                            end
                            
                        catch exception
                            
                            errorFileName =sprintf('%s.txt',wholeProFilePath);
                            % The error file name.
                            
                            msgString = getReport(exception,'extended','hyperlinks','off');
                            % Save error report as text file.
                            
                            fid = fopen(errorFileName,'w');
                            fprintf(fid, msgString);
                            fclose(fid);
                            
                            disp('Not processed');
                            % Save error message as a txt file
                        end
                        
                    elseif processFile == 0
                        disp('File skipped');
                        
                    end
                end
            end % frame.
        end
    end % Segment.
    
    %======================================================================
    % Save output .mat files.
    %======================================================================
    switch processingMode
        case 'concat'
            
            disp('...Saving file...')
            
            saveFileName = sprintf('%s/%spro/qualityControlData%s.mat',...
                iceBridgePath, dateStr(1:4),dateStr);
            
            save(saveFileName,...
                'dateNumVect',...
                'segmentNumVect',...
                'frameNumVect',...
                'traceNumVect',...
                'utcTimeVect',...
                'latitudeVect',...
                'longitudeVect',...
                'qcArray',...
                'atmDataAvailableVect',...
                'snowDataAvailableVect',...
                'minDistAtmPoint2NadirVect',...
                'minDistAtmPoint2NadirQcVect',...
                'maxDistAtmPoint2NadirVect',...
                'maxDistAtmPoint2NadirQcVect',...
                'meanAtmPitchVect',...
                'meanAtmPitchQcVect',...
                'meanAtmRollVect',...
                'meanAtmRollQcVect',...
                'hTopoVect',...
                'hTopoQcVect',...
                'interfaceOrderVect',...
                'interfaceOrderQcVect',...
                'minSnowDepthQcVect',...
                'maxSnowDepthQcVect',...
                'snrAirSnowVect',...
                'minSnrAirSnowQcVect',...
                'snrSnowIceVect',...
                'minSnrSnowIceQcVect',...
                'noiseAndClutterPowerVect',...
                'noisePowerThresholdVect',...
                'rangeBinAirSnowVect',...
                'rangeBinSnowIceVect',...
                'snowDepthVect',...
                'radarPrecisionVect',...
                'elevationMinVect',...
                'elevation05Vect',...
                'elevation50Vect',...
                'elevation95Vect',...
                'elevationMaxVect',...
                'elevationMeanVect',...
                'elevationStdVect',...
                'numAtmPointsVect',...
                'rangeBetweenInterfacesVect',...
                'rangeToAirSnowInterfaceVect',...
                'rangeToSnowIceInterfaceVect',...
                'fastTimeRangeVect',...
                'qualityFlagArray');
    end
    
    %======================================================================
    % Save qualFlagArray .mat files.
    disp('...Saving file...')
    %======================================================================
    switch processingMode
        case 'qualFlag'
            
            %--------------------------------------------------------------
            %% Save quality flag summary text file.
            %--------------------------------------------------------------
            qualityFlagArraySubset = qualityFlagArray(:,3:10);
            % The subset of the quality flag array.
            
            sumQfaRows = sum(qualityFlagArraySubset,1);
            % Sums the rows.
            
            sumQfaCols = sum(qualityFlagArraySubset,2);
            % Sums the cols.
            
            totalNumFiles = length(sumQfaCols);
            % The total number of files.
            
            numOfGoodFiles = sum(sumQfaCols==0);
            % The number of good files.
            
            numOfFlagFiles = sum(sumQfaCols>0);
            % The number of flagged files.
            
            % goodData          = '00000000';
            %
            % landOrIceberg     = '10000000';
            % unclassifiedError = '01000000';
            % lowSnr            = '00100000';
            % noGoodData        = '00010000';
            % missingData       = '00001000';
            % verticalStripes   = '00000100';
            % deconvolution     = '00000010';
            % coherentNoise     = '00000001';
            
            numLandOrIceberg     = sum(qualityFlagArraySubset(:,1));
            numUnclassifiedError = sum(qualityFlagArraySubset(:,2));
            numLowSnr            = sum(qualityFlagArraySubset(:,3));
            numNoGoodData        = sum(qualityFlagArraySubset(:,4));
            numMissingData       = sum(qualityFlagArraySubset(:,5));
            numVerticalStripes   = sum(qualityFlagArraySubset(:,6));
            numDeconvolution     = sum(qualityFlagArraySubset(:,7));
            numCoherentNoise     = sum(qualityFlagArraySubset(:,8));
            % The Number of files with the different quality flags.
            
            outputColumnNumber = [totalNumFiles;...
                numOfGoodFiles;...
                numOfFlagFiles;...
                numLandOrIceberg;...
                numUnclassifiedError;...
                numLowSnr;...
                numNoGoodData;...
                numMissingData;...
                numVerticalStripes;...
                numDeconvolution;...
                numCoherentNoise];
            % The output number array.
            
            outputColumnPercentages = 100*(outputColumnNumber./totalNumFiles);
            % The percentage array.
            
            comboOutArray = [outputColumnNumber,outputColumnPercentages]';
            % The combined array.
            
            %..............................................................
            %% Output text file.
            %..............................................................
            textFileName = sprintf('%s/%spro/qualityFlagSummary%s.txt',...
                iceBridgePath, dateStr(1:4),dateStr);
            
            fof1  = 'totalNumFiles (num,percent)           , %4.0f, %5.1f\n';
            fof2  = 'totalNumGoodFiles (num,percent)       , %4.0f, %5.1f\n';
            fof3  = 'totalNumFlaggedFiles (num,percent)    , %4.0f, %5.1f\n';
            fof4  = 'landOrIcebergFlagged (num,percent)    , %4.0f, %5.1f\n';
            fof5  = 'unclassifiedErrorFlagged (num,percent), %4.0f, %5.1f\n';
            fof6  = 'lowSnrFlagged (num,percent)           , %4.0f, %5.1f\n';
            fof7  = 'noGoodDataFlagged (num,percent)       , %4.0f, %5.1f\n';
            fof8  = 'missingDataFlagged (num,percent)      , %4.0f, %5.1f\n';
            fof9  = 'verticalStripesFlagged (num,percent)  , %4.0f, %5.1f\n';
            fof10 = 'deconvolutionFlagged (num,percent)    , %4.0f, %5.1f\n';
            fof11 = 'coherentNoiseFlagged (num,percent)    , %4.0f, %5.1f\n';
            % file format headers.
            
            fileOutputFormat =[fof1,fof2,fof3,fof4,fof5,fof6,fof7,fof8,fof9,fof10,fof11];
            % The file output format.
            
            fileID = fopen(textFileName,'w');
            fprintf(fileID, fileOutputFormat, comboOutArray(:));
            fclose(fileID);
            
            %--------------------------------------------------------------
            %% Save quality flag array.
            %--------------------------------------------------------------
            disp('...Saving file...')
            
            saveFileName = sprintf('%s/%spro/qualityFlagArray%s.mat',...
                iceBridgePath, dateStr(1:4),dateStr);
            
            save(saveFileName,'qualityFlagArray')
            %--------------------------------------------------------------
            
    end
    %======================================================================
    
end % date.
%==========================================================================


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% END OF CODE
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<