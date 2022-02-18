%%% Compile the DEM time series for the glacier of interest.

%INPUTS:
%1) Path containing all terrain parameter files.
%2) Two-letter abbreviation for study site. 
%3) Terrain parameter files with a consistent name format that starts with the date (YYYYMMDD) and
%includes the terrain parameter (DEM,aspect, slope, and ruggedness), and ends in .tif.

%OUTPUTS:
%1) matfile with all DEMs in the specified folder concatenated
%%%
%% Initialize
clearvars; close all;

DTM_path = '/Users/ellynenderlin/Research/NASA_CryoIdaho/glaciers/Wolverine/DEMs/';
abbrev = 'WG';
terrain = 'DEM';

%days of year
modays_norm = [31 28 31 30 31 30 31 31 30 31 30 31];
cumdays_norm = cumsum(modays_norm); cumdays_norm = [0 cumdays_norm(1:11)];
modays_leap = [31 29 31 30 31 30 31 31 30 31 30 31];
cumdays_leap = cumsum(modays_leap); cumdays_leap = [0 cumdays_leap(1:11)];


%% loop through all files for the specified terrain parameter and concatenate

%load the DEMs & convert dates to decimal years
cd_to_DEMs = ['cd ',DTM_path]; eval(cd_to_DEMs);
DEMs = dir(['*',terrain,'.tif']);
disp('Looping through terrain parameter files...');
for i = 1:length(DEMs)
    if contains(version,'2019')
        [Z(i).z,Z(i).R] = geotiffread(DEMs(i).name);
    else
        [Z(i).z,Z(i).R] = readgeoraster(DEMs(i).name);
    end
    
    %convert datestamp to decimal dates
    yr = str2double(DEMs(i).name(1:4));
    if mod(yr,4) == 0
        decidate = (cumdays_leap(str2double(DEMs(i).name(5:6)))+str2double(DEMs(i).name(7:8)))/sum(modays_leap);
    else
        decidate = (cumdays_norm(str2double(DEMs(i).name(5:6)))+str2double(DEMs(i).name(7:8)))/sum(modays_norm);
    end
    Z(i).deciyear = yr+decidate;
    
    %replace no data value with NaN
    Z(i).z(Z(i).z<=-9999) = NaN;
end

%create a data mask if working with the DEMs
if contains(terrain,'DEM')
    %identify the DEM with the smallest area
    for i = 1:length(Z)
        DEMarea(i) = Z(i).R.RasterExtentInWorldX.*Z(i).R.RasterExtentInWorldY;
    end
    maskref = find(DEMarea == min(DEMarea));
    
    %fill in gaps in the DEM
    DEMmask.z = ones(size(Z(maskref).z)); DEMmask.z(isnan(Z(maskref).z)) = 0;
    DEMmask.x = Z(maskref).R.XWorldLimits(1)+0.5*Z(maskref).R.CellExtentInWorldX:Z(maskref).R.CellExtentInWorldX:Z(maskref).R.XWorldLimits(end)-0.5*Z(maskref).R.CellExtentInWorldX; % get a vector of x coords
    DEMmask.y = Z(maskref).R.YWorldLimits(2)-0.5*Z(maskref).R.CellExtentInWorldY:-Z(maskref).R.CellExtentInWorldY:Z(maskref).R.YWorldLimits(1)+0.5*Z(maskref).R.CellExtentInWorldY; % get a vector of y coords
    clear maskref;
    %save
    save([abbrev,'-DEM-datamask.mat'],'DEMmask','-v7.3');
else
    load([abbrev,'-DEM-datamask.mat']);
end
%interpolate mask to each map & apply
[mask_xgrid,mask_ygrid] = meshgrid(DEMmask.x,DEMmask.y); 
for i = 1:length(Z)
    Z(i).x = Z(i).R.XWorldLimits(1)+0.5*Z(i).R.CellExtentInWorldX:Z(i).R.CellExtentInWorldX:Z(i).R.XWorldLimits(end)-0.5*Z(i).R.CellExtentInWorldX; % get a vector of x coords
    Z(i).y = Z(i).R.YWorldLimits(2)-0.5*Z(i).R.CellExtentInWorldY:-Z(i).R.CellExtentInWorldY:Z(i).R.YWorldLimits(1)+0.5*Z(i).R.CellExtentInWorldY;
    [terrain_xgrid,terrain_ygrid] = meshgrid(Z(i).x,Z(i).y); 
    interp_mask = interp2(mask_xgrid,mask_ygrid,DEMmask.z,terrain_xgrid,terrain_ygrid);
    interp_mask = round(interp_mask); interp_mask(isnan(interp_mask)) = 0; interp_mask = logical(interp_mask);
    Z(i).z(interp_mask==0) = NaN;
    clear terrain_*grid interp_mask;
end
%save the masked outputs
save([abbrev,'-',terrain,'-timeseries.mat'],'Z','-v7.3');
disp(['Concatenated ',terrain,' files']);
