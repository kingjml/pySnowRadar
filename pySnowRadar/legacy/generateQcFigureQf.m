function outputFigurePath = generateQcFigureQf(frame,surfaceAttributes,snow,snowQc,figFramePath)
% generateQcFigureQf - Generates quality control figures to assess the
% quality of the snow depth picks and ancillary information.
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    frame             - A structure containing the frame data.
%    surfaceAttributes - A structure containing the surfaceAttributes.
%    snow              - A structure containing the snow depth data.
%    snowQc            - A structure containing the snow qc.
%
% Outputs:
%    
%
% Example:
% -
%
% Other m-files required:
% - /nanmean.m
% - /stats/stats/quantile.m
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
% February 2016; Last revision: Base.

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% BEGIN CODE
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%frame
%surfaceAttributes
%snow
%snowQc


%==========================================================================
%% Setup.
%==========================================================================
close all

colordef white
% colordef black.

scrsz = get(0,'ScreenSize');
% Get the current screen dimensions.

%==========================================================================
%% Assign data from structures to temporary variables.
%==========================================================================

%--------------------------------------------------------------------------
% Frame variables.
%--------------------------------------------------------------------------
logEchogram = 20*log10(frame.Data);
% Extracts the echogram from the frame and converts it into log space.

logEchogram(~isfinite(logEchogram)) = nan;
% Changes non-finite numbers into nans.

[maxPowerSummed,maxPowerRangeBin] =  max(nanmean(logEchogram,2));
% Sums each trace in the echogram to find the range bin with the maximum
% power which will be the reference bin.

traceVect = 1:frame.echogramNumberOfCols;
% The vector for the number for traces.

rangeBinSeparation = 0.5*(frame.Time(2)-frame.Time(1))*snow.c;
% The fast time range bin separation.

%--------------------------------------------------------------------------
% SurfaceAttributes variables.
%--------------------------------------------------------------------------
meanAtmPitch = surfaceAttributes.meanAtmPitch;
% The mean ATM pitch for each snow radar measurement calculated using the
% ATM points closest to nadir (numAtmPointsForStats = 200).

meanAtmRoll  = surfaceAttributes.meanAtmRoll;
% The mean ATM roll for each snow radar measurement calculated using the
% ATM points closest to nadir (numAtmPointsForStats = 200).

h_topo       = surfaceAttributes.h_topo;
% The h_topo (95th - 5th elevation percentile) from each snow radar
% measurement calculated using the ATM points closest to nadir
% (numAtmPointsForStats = 200).

floeElevWav  = surfaceAttributes.floeElevWav;
% The wavelet floe elevation (multi-scale surface inflection point) from
% each snow radar measurement calculated using the ATM points closest to
% nadir (numAtmPointsForStats = 200).

h_topoMax = 0.5;
% Maximum h_topo = 0.5 m.

%..........................................................................
% Horizontal and vertical ATM distance.
%..........................................................................
numPercentiles = 17;
% The number for percentiles where the distance is sampled.

dtsAtmDistArray      = nan(numPercentiles,frame.echogramNumberOfCols);
% Setting up the distance from nadir array.

dtsAtmDistArrayColor = jet(numPercentiles);
% Setting up colorscale.

numTailPercentages = 4;
% The number of percentahes to sample the tails of the distribution.

dtsAtmElevArrayLower = nan(numTailPercentages,frame.echogramNumberOfCols);
% The lower percentiles: [0, 1, 2.5, 5, 10]

dtsAtmElevArrayUpper = nan(numTailPercentages,frame.echogramNumberOfCols);
% The upper percentiles: [90, 95, 97.5, 99, 100]

for ind = traceVect
    if isstruct(surfaceAttributes.dtsAtmDist{1,ind})
        
        dtsAtmDistArray(1 ,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile0000;
        dtsAtmDistArray(2 ,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile0010;
        dtsAtmDistArray(3 ,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile0025;
        dtsAtmDistArray(4 ,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile0050;
        dtsAtmDistArray(5 ,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile0100;
        dtsAtmDistArray(6 ,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile0200;
        dtsAtmDistArray(7 ,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile0300;
        dtsAtmDistArray(8 ,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile0400;
        dtsAtmDistArray(9 ,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile0500;
        dtsAtmDistArray(10,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile0600;
        dtsAtmDistArray(11,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile0700;
        dtsAtmDistArray(12,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile0800;
        dtsAtmDistArray(13,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile0900;
        dtsAtmDistArray(14,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile0950;
        dtsAtmDistArray(15,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile0975;
        dtsAtmDistArray(16,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile0990;
        dtsAtmDistArray(17,ind) = surfaceAttributes.dtsAtmDist{1,ind}.Percentile1000;
        
        dtsAtmElevArrayLower(1,ind) = surfaceAttributes.dtsAtmElev{1,ind}.Percentile0010;
        dtsAtmElevArrayLower(2,ind) = surfaceAttributes.dtsAtmElev{1,ind}.Percentile0025;
        dtsAtmElevArrayLower(3,ind) = surfaceAttributes.dtsAtmElev{1,ind}.Percentile0050;
        dtsAtmElevArrayLower(4,ind) = surfaceAttributes.dtsAtmElev{1,ind}.Percentile0100;
        
        dtsAtmElevArrayUpper(1,ind) = surfaceAttributes.dtsAtmElev{1,ind}.Percentile0900;
        dtsAtmElevArrayUpper(2,ind) = surfaceAttributes.dtsAtmElev{1,ind}.Percentile0950;
        dtsAtmElevArrayUpper(3,ind) = surfaceAttributes.dtsAtmElev{1,ind}.Percentile0975;
        dtsAtmElevArrayUpper(4,ind) = surfaceAttributes.dtsAtmElev{1,ind}.Percentile0990;
          
    end
end
% This array contains the percentiles of the ATM points closest to nadir
% (numAtmPointsForStats = 200)
%..........................................................................

%--------------------------------------------------------------------------
% Snow variables.
%--------------------------------------------------------------------------
rangeBinVectAirSnow = snow.rangeBinVectAirSnow;
% The air/snow range bin location.

rangeBinVectSnowIce = snow.rangeBinVectSnowIce;
% The snow/ice range bin location.

snowDepth = snow.depth;
% The distance between first and maximum return converted into a snow
% depth.
%--------------------------------------------------------------------------



%==========================================================================
% Figure 1: Echogram.
ax(1) = subplot(10,2,1:8);
set(gcf,'Visible','off')
%==========================================================================
imagesc(traceVect,[],logEchogram)
% Plots the echogram.

hold on

%--------------------------------------------------------------------------
% Plot the first return.
%--------------------------------------------------------------------------
rangeBinVectAirSnowDisp = rangeBinVectAirSnow;
rangeBinVectAirSnowDisp(h_topo>h_topoMax) = nan;
% Only display the picks that meet the h_topo condition.

plot(traceVect,rangeBinVectAirSnowDisp,'-',...
    'color','r',...
    'LineWidth',1)

%--------------------------------------------------------------------------
% Plot the maximum return.
%--------------------------------------------------------------------------
rangeBinVectSnowIceDisp = rangeBinVectSnowIce;
rangeBinVectSnowIceDisp(h_topo>h_topoMax) = nan;
% Only display the picks that meet the h_topo condition.

plot(traceVect,rangeBinVectSnowIceDisp,'-',...
    'color','k',...
    'LineWidth',1)

%--------------------------------------------------------------------------
% Plot the color axis.
%--------------------------------------------------------------------------
caxisMin = quantile(quantile(logEchogram,0.50),0.50);
% The caxis min.

caxisMax = quantile(max(logEchogram),0.95);
% The caxis max.

if caxisMin<caxisMax
    caxis([caxisMin caxisMax])
end
% set caxis.

colormap(parula(1000))
% Set the colormap to the parula colormap.

%--------------------------------------------------------------------------
% Set x and y axis limits.
%--------------------------------------------------------------------------
xlim([traceVect(1) traceVect(end)])
% The x limits.

deltaRangeMetresUpper = 2;
% The number of metres above the maxPowerRangeBin to crop the ylimit to.

deltaRangeBinsUpper = ceil(deltaRangeMetresUpper/rangeBinSeparation);
% The number of range bins above the maxPowerRangeBin to crop the ylimit to.

deltaRangeMetresLower = 1;
% The number of metres below the maxPowerRangeBin to crop the ylimit to.

deltaRangeBinsLower = ceil(deltaRangeMetresLower/rangeBinSeparation);
% The number of range bins below the maxPowerRangeBin to crop the ylimit to.

if frame.echogramNumberOfRows<=2000
    ylim([maxPowerRangeBin-deltaRangeBinsUpper,maxPowerRangeBin+deltaRangeBinsLower])
    % The y limits.
    
elseif frame.echogramNumberOfRows > 2000
    ylim([1,frame.echogramNumberOfRows])
    % The y limits.
end

grid on
box on

%--------------------------------------------------------------------------
% Titles and labels.
%--------------------------------------------------------------------------
saveFileName = sprintf('%s_%s_%02i_%03i',frame.radarType,frame.dateStr,frame.segmentNumber,frame.frameNumber);
% generate the save file name and plot title.

titleStr = sprintf('%s - Qf: %s',saveFileName,frame.qualityFlag);
% The title string.

title(titleStr,...
    'fontSize',12,...
    'BackgroundColor',[1 1 1],...
    'Interpreter', 'none');
% Set the title with a white background.

ylabel('Range bin','fontSize',10)
% Sets the y label.

set(gca,'XTickLabel',[])
% Turns the x label off.
%--------------------------------------------------------------------------


%==========================================================================
%% Figure 2: The snow depth.
ax(2) = subplot(10,2,9:20);
%==========================================================================
imagesc(traceVect,[],logEchogram)
% Plots the echogram.

if caxisMin<caxisMax
    caxis([caxisMin caxisMax])
end

ylabel('Range bin','fontSize',10)
% Sets the y label.

xlabel('Trace number','fontSize',10)
% Sets the x label.?

%==========================================================================
% Output figure.
%==========================================================================
outputFigurePath = sprintf('%s/%sQf.png',figFramePath,saveFileName);

screenRatio = 1280/800;
% the screen ratio.

pageLength = 11;
pageWidth  = pageLength/screenRatio ;
% The page dimensions.

set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 pageLength pageWidth]);
print(gcf,outputFigurePath,'-r300','-dpng');


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% END OF CODE
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end