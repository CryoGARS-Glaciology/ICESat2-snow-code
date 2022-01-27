function [ztruth,differences,rmsez] = extract_icesat2_vertical_errors(icesat2, elevations, R2)
% Function extract_icesat2_vertical_errors computes the vertical offsets
% beteen horizontally-coregistered ICESat-2 tracks and a reference DTM
% INPUTS: icesat2 = a csv file with icesat 2 elevations created using the
%                       h5 to csv jupyter notebook
%              elevations = the matrix created using readgeoraster()
%              R2 = the cell map refernce created as the second output in
%                       readgeoraster()
% OUTPUTS:      ztruth = reference elevation transect
%               differences = difference in elevation between the icesat-2
%                       footprints and their corresponding DTM locations
%               rmsez = the root mean squared elevation difference between the 
%                       icesat-2 footprints and corresponding DTM locations

% Created 19 October 2020 by Colten Elkin (coltenelkin@u.boisestate.edu)
% last modified 08 June 2021 by Ellyn Enderlin (ellynenderlin@boisestate.edu)

%filter reference DTM elevations
% elevations(elevations < -10) = nan; 
elevations(elevations > 10000) = nan; % throw out bad data (clouds)

%read ICESat-2 transect
T = readtable(icesat2);
%if you have vertically coregistered the ICESat2 data, use the adjusted elevations
if ismember('Elevation_Coregistered', T.Properties.VariableNames)
    zmod = T.Elevation_Coregistered(:);
else
    zmod = T.Elevation(:); % save the 'model' elevations (icesat-2 elevations)
end
easts = T.Easting(:); % pull out the easting values
norths = T.Northing(:); % pull out the northings
footwidth = 11; % approx. width of icesat2 shot footprint in meters

% initialize empty matrices
theta = zeros(size(norths)); 
xs = {};
ys = {};
xpoly = nan([1,5]);
ypoly = nan([1,5]);

%create polygons of ICESat-2 footprints
for r = 1:length(theta)
    if r == length(theta)
        theta(r) = abs(atan((norths(r) - norths(r-1))/(easts(r) - easts(r-1)))); % trig to get angle theta along-track
    else
        theta(r) = abs(atan((norths(r+1) - norths(r))/(easts(r+1) - easts(r)))); % trig to get angle theta along-track
    end
    
    if contains(icesat2, 'ATL08') % ATL08 commands   
        % get the x and y vectors to form the polygon
        xpoly(1) = easts(r) + (footwidth/2) - footwidth/2*cos((pi/2) - theta(r)); % calculate the 4 corners in the x direction
        xpoly(2) = easts(r) + (footwidth/2) + footwidth/2*cos((pi/2) - theta(r));
        xpoly(3) = easts(r) - (footwidth/2) + footwidth/2*cos((pi/2) - theta(r));
        xpoly(4) = easts(r) - (footwidth/2) - footwidth/2*cos((pi/2) - theta(r));
        xs{r} = [xpoly(1), xpoly(2), xpoly(3), xpoly(4), xpoly(1)]; % save the corners as a vector in the x-es cell array
        
        ypoly(1) = norths(r) - 50 - footwidth/2*cos((pi/2) - theta(r)); % calculate the 4 corners in the y direction
        ypoly(2) = norths(r) - 50 + footwidth/2*cos((pi/2) - theta(r));
        ypoly(3) = norths(r) + 50 + footwidth/2*cos((pi/2) - theta(r));
        ypoly(4) = norths(r) + 50 - footwidth/2*cos((pi/2) - theta(r));
        ys{r} = [ypoly(1), ypoly(2), ypoly(3), ypoly(4), ypoly(1)]; % save the corners as a vector in the y-s cell array
        
    elseif contains(icesat2, 'ATL06') % ATL06 commands
        % get the x and y vectors to form the polygon
        xpoly(1) = easts(r) + (footwidth/2) - footwidth/2*cos((pi/2) - theta(r)); % calculate the 4 corners in the x direction
        xpoly(2) = easts(r) + (footwidth/2) + footwidth/2*cos((pi/2) - theta(r));
        xpoly(3) = easts(r) - (footwidth/2) + footwidth/2*cos((pi/2) - theta(r));
        xpoly(4) = easts(r) - (footwidth/2) - footwidth/2*cos((pi/2) - theta(r));
        xs{r} = [xpoly(1), xpoly(2), xpoly(3), xpoly(4), xpoly(1)]; % save the corners as a vector in the x-es cell array
        
        ypoly(1) = norths(r) - 20 - footwidth/2*cos((pi/2) - theta(r)); % calculate the 4 corners in the y direction
        ypoly(2) = norths(r) - 20 + footwidth/2*cos((pi/2) - theta(r));
        ypoly(3) = norths(r) + 20 + footwidth/2*cos((pi/2) - theta(r));
        ypoly(4) = norths(r) + 20 - footwidth/2*cos((pi/2) - theta(r));
        ys{r} = [ypoly(1), ypoly(2), ypoly(3), ypoly(4), ypoly(1)]; % save the corners as a vector in the y-s cell array
    end
end

%define the reference elevation data
x = R2.XWorldLimits(1)+0.5*R2.CellExtentInWorldX:R2.CellExtentInWorldX:R2.XWorldLimits(2)-0.5*R2.CellExtentInWorldX; % get a vector of x coords
if strcmp(R2.ColumnsStartFrom,'north')
    y = R2.YWorldLimits(2)-0.5*R2.CellExtentInWorldY:-R2.CellExtentInWorldY:R2.YWorldLimits(1)+0.5*R2.CellExtentInWorldY; 
else
    y = R2.YWorldLimits(1)+0.5*R2.CellExtentInWorldY:R2.CellExtentInWorldY:R2.YWorldLimits(2)-0.5*R2.CellExtentInWorldY; 
end
[xgrid, ygrid] = meshgrid(x, y); % create grids of each of the x and y coords
elevation_report = zeros([1, length(xs)]);

%identify the reference elevation points in each ICESat2 footprint
for t = 1:length(xs)
    xv = xs{t}; % bounding box x vector
    yv = ys{t}; % bounding box y vector
    
    % first trimming
    in = inpolygon(xgrid, ygrid, xv, yv); % get logical array of in values
    pointsinx = xgrid(in); % save x locations
    pointsiny = ygrid(in); % save y locations
    elevationsin = elevations(in); % save elevations
    elevation_report(t) = nanmean(elevationsin); %NEED TO MODIFY... 
    %%FIGURE OUT HOW TO WEIGHT ACROSS-TRACK & HOW TO FIT A SLOPE LIKE
    %%h_te_best_fit VARIABLE SAVED WITH ATL08
end
ztruth = elevation_report(:);

%calculate the elevation residuals
differences = zmod - ztruth;
% differences(differences > 80) = nan; differences(differences < -80) = nan;
rmsez = sqrt(nansum((differences).^2)./length(differences));
