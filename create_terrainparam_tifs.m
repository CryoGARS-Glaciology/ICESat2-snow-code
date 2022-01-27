%%% Compute slope and aspect from DEMs.

%INPUTS:
%1) Path containing all DEMs.
%2) UTM zone number and EPSG code.

%OUTPUTS:
%1) tifs for elevation, slope, and aspect in the native (UTM) coordinates

%DEPENDENCIES:
%1) utm2ll to transform UTM coordinates to latitude and longitude in order
%to use Matlab's built-in gradientm function. Path specified in code.
%%%
%% Initialize
clear all; close all;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/');

%specify the directory and UTM zone
cd /users/ellynenderlin/Research/NASA_CryoIdaho/glaciers/Wolverine/DEMs/
UTMzone = 6; EPSGcode = 32606;

%% loop through the DEMs, calculating slope and aspect & resaving as geotiffs
DEMs = dir('*.tif');
for j = 1:length(DEMs)
    disp(DEMs(j).name);
    [Z,R] = geotiffread(DEMs(j).name); Z(Z<0) = NaN;
    [xgrid,ygrid] = meshgrid(R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2),R.YWorldLimits(2):-R.CellExtentInWorldY:R.YWorldLimits(1));
    
    [LAT,LON]=utm2ll(xgrid,ygrid,UTMzone);
    lat_spacing = nanmean(diff(LAT(:,1),1,1)); lon_spacing = nanmean(diff(LON(1,:),1,2)); 
    latlim = [min(min(LAT)) max(max(LAT))]; lonlim = [min(min(LON)) max(max(LON))];
    R2 = georefcells(latlim,lonlim,R.RasterSize);
    [ASPECT,SLOPE,~,~] = gradientm(double(Z),R2); %aspect = degrees clockwise from north, slope = degrees up from horizontal
    ASPECT(isnan(Z)) = NaN; SLOPE(isnan(Z)) = NaN;
    
    %plot the data
    figure('position',[50 50 1200 400]);
    subplot(1,3,1); %elevation
    imagesc(R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2),...
        R.YWorldLimits(2):-R.CellExtentInWorldY:R.YWorldLimits(1),Z); 
    axis xy equal; colormap parula; colorbar; 
    subplot(1,3,2); %slope
    imagesc(R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2),...
        R.YWorldLimits(2):-R.CellExtentInWorldY:R.YWorldLimits(1),SLOPE); 
    axis xy equal; colormap parula; 
    set(gca,'clim',[0 30]); colorbar; 
    subplot(1,3,3); %aspect
    imagesc(R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2),...
        R.YWorldLimits(2):-R.CellExtentInWorldY:R.YWorldLimits(1),ASPECT); 
    axis xy equal; colormap parula; colorbar; 
    drawnow;
    
    %save the data
    geotiffwrite([DEMs(j).name(1:end-4),'2.tif'],Z,R,'CoordRefSysCode',['EPSG:',num2str(EPSGcode),'']); %elevations
    geotiffwrite([DEMs(j).name(1:end-4),'-slope.tif'],SLOPE,R,'CoordRefSysCode',['EPSG:',num2str(EPSGcode),'']); %slope
    geotiffwrite([DEMs(j).name(1:end-4),'-aspect.tif'],ASPECT,R,'CoordRefSysCode',['EPSG:',num2str(EPSGcode),'']); %aspect
    disp('... terrain parameters saved');
    clear R* Z SLOPE ASPECT *grid LAT LON lat* lon*;
end