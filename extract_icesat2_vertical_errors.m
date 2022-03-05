function [ztruth,differences,rmsez] = extract_icesat2_vertical_errors(icesat2, elevations, R2)
% Function extract_icesat2_vertical_errors computes the vertical offsets
% beteen horizontally-coregistered ICESat-2 tracks and a reference DTM
% INPUTS: icesat2 = horizontally-coregistered csv file with icesat 2 elevations
%              elevations = the matrix created using readgeoraster()
%              R2 = the cell map refernce created as the second output in
%                       readgeoraster()
% OUTPUTS:      ztruth = reference elevation transect
%               differences = difference in elevation between the icesat-2
%                       footprints and their corresponding DTM locations
%               rmsez = the root mean squared elevation difference between the 
%                       icesat-2 footprints and corresponding DTM locations

% last modified 03 March 2022 by Ellyn Enderlin (ellynenderlin@boisestate.edu)

%specify ICESat-2 footprint width & length
footwidth = 11; % approx. width of icesat2 shot footprint in meters
if contains(icesat2, 'ATL08') % ATL08 commands
    default_length = 100; % approx. length of icesat2 shot footprint in meters
elseif contains(icesat2, 'ATL06') % ATL06 commands
    default_length = 40; % approx. length of icesat2 shot footprint in meters
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

% initialize empty matrices
theta = NaN(length(norths),2); 

%create polygons of ICESat-2 footprints
for r = 1:length(theta)
    %calculate the footprint orientation
    if r == 1 %only angle pointed forwards
        theta(r,1) = atan2d((norths(r+1)-norths(r)),(easts(r+1)-easts(r))); theta(r,2) = theta(r,1);
        footlength = default_length;
    elseif r == length(theta) %only angle pointed backwards
        theta(r,1) = atan2d((norths(r)-norths(r-1)),(easts(r)-easts(r-1))); theta(r,2) = theta(r,1);
        footlength = default_length;
    else %calculate angles in each direction
        theta(r,1) = atan2d((norths(r)-norths(r-1)),(easts(r)-easts(r-1)));
        theta(r,2) = atan2d((norths(r+1)-norths(r)),(easts(r+1)-easts(r)));
        footlength = sqrt((norths(r+1)-norths(r-1)).^2 + (easts(r+1)-easts(r-1)).^2)/2;
    end
    
    %find box edges along the RGT
    back_x = easts(r)-(footlength/2)*cosd(theta(r,1)); back_y = norths(r)-(footlength/2)*sind(theta(r,1));
    front_x = easts(r)+(footlength/2)*cosd(theta(r,2)); front_y = norths(r)+(footlength/2)*sind(theta(r,2));
    
    %find box edges perpendicular to the centroid
    xc(r,1) = easts(r)+(footwidth/2)*cosd(nanmean(theta(r,:))+90); yc(r,1) = norths(r)+(footwidth/2)*sind(nanmean(theta(r,:))+90);
    xc(r,4) = easts(r)+(footwidth/2)*cosd(nanmean(theta(r,:))-90); yc(r,4) = norths(r)+(footwidth/2)*sind(nanmean(theta(r,:))-90);
    
    %solve for corner coordinates
    xc(r,2) = back_x+(footwidth/2)*cosd(theta(r,1)+90); yc(r,2) = back_y+(footwidth/2)*sind(theta(r,1)+90);
    xc(r,3) = back_x+(footwidth/2)*cosd(theta(r,1)-90); yc(r,3) = back_y+(footwidth/2)*sind(theta(r,1)-90);
    xc(r,5) = front_x+(footwidth/2)*cosd(theta(r,2)-90); yc(r,5) = front_y+(footwidth/2)*sind(theta(r,2)-90);
    xc(r,6) = front_x+(footwidth/2)*cosd(theta(r,2)+90); yc(r,6) = front_y+(footwidth/2)*sind(theta(r,2)+90);
    clear back_* front_*;
end

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
for k = 1:length(xc)
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
