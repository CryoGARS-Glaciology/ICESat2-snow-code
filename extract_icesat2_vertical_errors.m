function [ztruth,differences,rmsez] = extract_icesat2_vertical_errors(icesat2, elevations, R2)
% Function EXTRACT_ICESAT2_VERTICAL_ERRORS computes the vertical offsets
% beteen horizontally-coregistered ICESat-2 tracks and relevation elevation
% dataset(s).
% INPUTS: icesat2 = csv file(s) with ICESat-2 data saved as column vectors
%      elevations = the reference elevation matrix 
%              R2 = the cell map reference for the reference DTM
% OUTPUTS:      ztruth = reference elevation transect
%               differences = difference in elevation between the icesat-2
%                       segments and their corresponding DTM locations
%               rmsez = the root mean squared elevation difference between the 
%                       icesat-2 segments and corresponding DTM locations

% last modified 03 March 2022 by Ellyn Enderlin (ellynenderlin@boisestate.edu)

%specify ICESat-2 footprint width & length
footwidth = 11; % approx. width of icesat2 shot footprint in meters
if contains(icesat2, 'ATL08') % ATL08 commands
    default_length = 100; % approx. length of icesat2 shot footprint in meters
    ATL0X = 8; %dataset flag for footprint delineation
elseif contains(icesat2, 'ATL06') % ATL06 commands
    default_length = 40; % approx. length of icesat2 shot footprint in meters
    ATL0X = 6; %dataset flag for footprint delineation
end

%filter reference DTM elevations
% elevations(elevations < -10) = nan; 
elevations(elevations > 10000) = nan; % throw out bad data (clouds)

%if you have vertically coregistered the ICESat2 data, use the adjusted elevations
T = readtable(icesat2);
if ismember('Elevation_Coregistered', T.Properties.VariableNames)
    zmod = T.Elevation_Coregistered(:);
else
    zmod = T.Elevation(:); % save the 'model' elevations (icesat-2 elevations)
end
easts = T.Easting(:); % pull out the easting values
norths = T.Northing(:); % pull out the northings
end_flag = [1; zeros(length(zmod)-2,1); 1];

%create polygons of ICESat-2 segments
[xc,yc,~] = ICESat2_FootprintCorners(norths,easts,ATL0X,end_flag);

%define the reference elevation data
x = R2.XWorldLimits(1)+0.5*R2.CellExtentInWorldX:R2.CellExtentInWorldX:R2.XWorldLimits(2)-0.5*R2.CellExtentInWorldX; % get a vector of x coords
if strcmp(R2.ColumnsStartFrom,'north')
    y = R2.YWorldLimits(2)-0.5*R2.CellExtentInWorldY:-R2.CellExtentInWorldY:R2.YWorldLimits(1)+0.5*R2.CellExtentInWorldY; 
else
    y = R2.YWorldLimits(1)+0.5*R2.CellExtentInWorldY:R2.CellExtentInWorldY:R2.YWorldLimits(2)-0.5*R2.CellExtentInWorldY; 
end
[xgrid, ygrid] = meshgrid(x, y); % create grids of each of the x and y coords
elevation_report = zeros([1, size(xc,1)]);

%identify the reference elevation points in each ICESat2 footprint
for k = 1:size(xc,1)
    xv = [xc(k,:) xc(k,1)]; % bounding box x vector
    yv = [yc(k,:) yc(k,1)]; % bounding box y vector
    
    %data in the footprint
    in = inpolygon(xgrid, ygrid, xv, yv); % get logical array of in values
    elevationsin = elevations(in); % isolate elevations within the footprint
    elevation_report(k) = nanmean(elevationsin); %export the mean
    clear xv yv in elevationsin;
end
ztruth = elevation_report(:);

%calculate the elevation residuals
differences = zmod - ztruth;
% differences(differences > 80) = nan; differences(differences < -80) = nan;
rmsez = sqrt(nansum((differences).^2)./length(differences));
