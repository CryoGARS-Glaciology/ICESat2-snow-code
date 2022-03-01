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

%specify the file path and naming conventions
DTM_path = '/Users/ellynenderlin/Research/NASA_CryoIdaho/glaciers/Wolverine/DEMs/';
abbrev = 'WG'; %site abbreviation
terrain = 'DEM'; %terrain parameter to concatenate into a .mat file

%For glacierized setting, specify the site-specific glacier & snow info:
%glacier outline needed to exclude all on-ice locations from coreg
%snowline to make sure small unidentified ice fields are not included
S = shaperead('/Users/ellynenderlin/Research/NASA_CryoIdaho/glaciers/Wolverine/ROIs/Wolverine-2018-outline-UTM06N.shp'); %glacier outline in UTM coords (same as DEMs)
snowline = 1235; %Wolverine Glacier = 1235; 

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

%determine lowest-resolution DEM for coregistration purposes
if contains(terrain,'DEM')
    for i = 1:length(Z)
        pixelsz(i) = abs(Z(i).x(1)-Z(i).x(2))*abs(Z(i).y(1)-Z(i).y(2));
    end
    refDEMind = find(pixelsz==max(pixelsz),1,'first');
    [refDEMx,refDEMy] = meshgrid(Z(refDEMind).x,Z(refDEMind).y);
    glacier = inpolygon(refDEMx,refDEMy,S.X,S.Y);
end

%create a data mask if working with the DEMs
maskfile = dir(['*-DEM-datamask.mat']);
if contains(terrain,'DEM') && isempty(maskfile)
    disp('creating a mask of the minimum DEM extent');
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
    disp('reloading existing mask');
    load([abbrev,'-DEM-datamask.mat']);
end

%interpolate mask to each map & apply
[mask_xgrid,mask_ygrid] = meshgrid(DEMmask.x,DEMmask.y); 
for i = 1:length(Z)
    %apply mask
    Z(i).x = Z(i).R.XWorldLimits(1)+0.5*Z(i).R.CellExtentInWorldX:Z(i).R.CellExtentInWorldX:Z(i).R.XWorldLimits(end)-0.5*Z(i).R.CellExtentInWorldX; % get a vector of x coords
    Z(i).y = Z(i).R.YWorldLimits(2)-0.5*Z(i).R.CellExtentInWorldY:-Z(i).R.CellExtentInWorldY:Z(i).R.YWorldLimits(1)+0.5*Z(i).R.CellExtentInWorldY;
    [terrain_xgrid,terrain_ygrid] = meshgrid(Z(i).x,Z(i).y); 
    interp_mask = interp2(mask_xgrid,mask_ygrid,DEMmask.z,terrain_xgrid,terrain_ygrid);
    interp_mask = round(interp_mask); interp_mask(isnan(interp_mask)) = 0; interp_mask = logical(interp_mask);
    Z(i).z(interp_mask==0) = NaN;
    
    %vertically adjust so all DEMs are decently coregistered
    if contains(terrain,'DEM')
        if i ~= refDEMind
            [coregDEMx,coregDEMy] = meshgrid(Z(i).x,Z(i).y);
            DEMinterp = interp2(coregDEMx,coregDEMy,double(Z(i).z),refDEMx,refDEMy);
            DEMdiff = DEMinterp - Z(refDEMind).z;
            coregz = nanmedian(DEMdiff(Z(refDEMind).z<snowline & glacier==0));
            Z(i).z = Z(i).z - coregz; Z(i).z_coreg = DEMdiff;
            clear coreg* DEMinterp DEMdiff;
        end
    end
    clear terrain_*grid interp_mask;
end
%save the masked outputs
save([abbrev,'-',terrain,'-timeseries.mat'],'Z','-v7.3');
disp(['Concatenated ',terrain,' files']);
