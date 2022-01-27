%calculate the EGM96 geoid height for RCEW
clearvars; close all;
addpath('/Users/ellynenderlin/mfiles/general/');

% %Reynolds Creek Experimental Watershed, Idaho, from the Idaho Lidar
% %Consortium (NAVD88 vertical reference, NAD83 UTM11N horizontal reference)
% DTM_path = '/Users/ellynenderlin/Research/NASA_CryoIdaho/mountains/RCEW/DEMs/';
%Dry Creek Experimental Watershed, Idaho, from the Idaho Lidar
%Consortium (NAVD88 vertical reference, NAD83 UTM11N horizontal reference)
DTM_path = '/Users/ellynenderlin/Research/NASA_CryoIdaho/mountains/DC/DEMs/';

%site abbreviation for file names
site_abbrev = 'DC';

%DTM name
% DTM_name = '2014_Lidar_Derived_1m_DEM_RCEW.tif'; %Reynolds Creek 
DTM_name = '2007_DryCreekDEM-NAVD88-NAD83.tif'; %Dry Creek
DTM_outputname = '2007_DryCreekDEM-ellipsoidalNAD83.tif';

%% extract cartesian coordinates for use in NOAA NGS GEOID12b INTG.exe package

inc = 20; %sparse grid increment

%load the DTM
cd_to_DTM = ['cd ',DTM_path]; eval(cd_to_DTM);
[myDTM,myRef] = geotiffread(DTM_name);
x = myRef.XWorldLimits(1):inc*myRef.CellExtentInWorldX:myRef.XWorldLimits(2);
y = myRef.YWorldLimits(2):-inc*myRef.CellExtentInWorldY:myRef.YWorldLimits(1);
[xgrid,ygrid] = meshgrid(x,y);

%convert UTM 11N coordinates to Cartesian coordinates
coords = [];
counter = 0; %tic;
for i = 1:size(xgrid,1)
    counter = counter+1;
    for j = 1:size(xgrid,2)
        [lat,lon]=utm2ll(xgrid(i,j),ygrid(i,j),11,'nad83'); %set-up for UTM11N in NAD83! (Idaho Lidar Consortium format)
        coords = [coords; lat 360+lon];
        clear lat lon;
    end
    
    %provide a time update
    if counter == 1000
        counter = 0; 
%         toc; tic;
    end
end
% save RCEW-DTM-coords.txt coords -ascii;
% dlmwrite('RCEW-DTM-coords.txt',coords,'precision','%.4f','delimiter','\t');
fnm = fullfile(DTM_path,[site_abbrev,'-DTM-coords.txt']);
fid = fopen(fnm,'wt');
fprintf(fid,'%.6f %.6f\n',coords');
fclose(fid);
clear DTM_name;
disp('coordinates saved as an ascii text file');
disp('now go to PC & plug coords into https://www.ngs.noaa.gov/GEOID/GEOID12B/GEOID12B_CONUS.shtml INTG.exe');


%% interpolate sparse geoid12b heights to full DTM grid
coords = readmatrix([site_abbrev,'-DTM-coords.txt']);
linref = 1:length(coords); %create a vector reference for the extracted coordinates from the full grid

%load the geoid heights from the INTG.exe output csv
geoid_output = readmatrix([site_abbrev,'-NAVD88-heights.txt']);
geoid_vector = geoid_output(:,7);

%extract the full coordinate grid
xfull = myRef.XWorldLimits(1)+0.5*myRef.CellExtentInWorldX:myRef.CellExtentInWorldX:myRef.XWorldLimits(2)-0.5*myRef.CellExtentInWorldX;
yfull = myRef.YWorldLimits(2)-0.5*myRef.CellExtentInWorldY:-myRef.CellExtentInWorldY:myRef.YWorldLimits(1)+0.5*myRef.CellExtentInWorldY;
[xfull_grid,yfull_grid] = meshgrid(xfull,yfull);

%convert ascii text geoid heights to a matrix
[row,col] = ind2sub(size(xgrid),linref); %linear indices to subscripts
for i = 1:length(linref)
    geoid_matrix(row(i),col(i)) = geoid_vector(linref(i));
end

%interpolate from the sparse grid to the full coordinate matrix
geoidfull = interp2(xgrid,ygrid,geoid_matrix,xfull_grid,yfull_grid,'markima');
geoidsmooth = imgaussfilt(geoidfull,100); %smooth strange step-like features

%adjust to ellipsoidal elevations by adding the geoid to the DTM
Z = myDTM+geoidsmooth;
geotiffwrite(DTM_outputname, Z, myRef, 'CoordRefSysCode', 'EPSG:26911'); %specified for NAD83 UTM11N
disp('vertically adjusted to ellipsoidal heights');

%% compute EGM96 heights
% DTM_name = '2014_Lidar_Derived_1m_DEM_RCEW-WGS84(OLD).tif';
% [myDTM,myRef] = geotiffread(DTM_name);
% lon = myRef.LongitudeLimits(1)+0.5*myRef.CellExtentInLongitude:myRef.CellExtentInLongitude:myRef.LongitudeLimits(2)-0.5*myRef.CellExtentInLongitude;
% lat = myRef.LatitudeLimits(2)-0.5*myRef.CellExtentInLatitude:-myRef.CellExtentInLatitude:myRef.LatitudeLimits(1)+0.5*myRef.CellExtentInLatitude;
% [lon_grid,lat_grid] = meshgrid(lon,lat);
% EGM96height = geoidheight(lat_grid,lon_grid);
% 
% save('RCEW_2014_Lidar_Derived_1m_DEM-geoidheights.mat','EGM96height');
% % geotiffwrite('RCEW_2014_Lidar_Derived_1m_DEM-geoidheights.tif', EGM96height, myRef, 'CoordRefSysCode', 'EPSG:4326');
% disp('EGM96 geiod heights saved');

