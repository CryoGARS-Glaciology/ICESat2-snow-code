function [params, range, SDev] = calc_icesat2_params(icesat2, tif, R2)
% Function COREGISTER_ICESAT2 coregisters icesat-2 data with a corresponding digital
% terrain model 
% INPUTS: icesat2 = a csv file with icesat 2 elevations created using the
%                       h5 to csv jupyter notebook
%             tif = the reference matrix or structure of the terrain
%                       parameter (slope, aspect, etc)
%              R2 = the cell map reference
% OUTPUTS: params = the terrain parameter from the tif averaged over the
%                   bounding box of the icesat2 footprint. Reported as a
%                   column vector matching the icesat2 input points
%           range = the range (max - min) of the parameter in the footprint
%            SDev = the standard deviation of the parameter in the
%                   footprint

% Created 30 November 2020 by Colten Elkin (coltenelkin@u.boisestate.edu)
% Last modified 09 June 2021 by Ellyn Enderlin (ellynenderlin@boisestate.edu)

T = readtable(icesat2);
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

%extract data from the input geotiff if over static terrain (ATL08) or from
%the geotiff timeseries if over a temporally-evolving surface (ATL06)
if nargin == 3
    %define the reference elevation data
    x = R2.XWorldLimits(1)+0.5*R2.CellExtentInWorldX:R2.CellExtentInWorldX:R2.XWorldLimits(end)-0.5*R2.CellExtentInWorldX; % get a vector of x coords
    y = R2.YWorldLimits(2)-0.5*R2.CellExtentInWorldY:-R2.CellExtentInWorldY:R2.YWorldLimits(1)+0.5*R2.CellExtentInWorldY; % get a vector of y coords
    [xgrid, ygrid] = meshgrid(x, y); % create grids of each of the x and y coords
    parameter_report = NaN([1, length(xs)]);
    
    %identify the reference elevation points in each ICESat2 footprint
    for t = 1:length(xs)
        xv = xs{t}; % bounding box x vector
        yv = ys{t}; % bounding box y vector
        
        % first trimming
        in = inpolygon(xgrid, ygrid, xv, yv); % get logical array of in values
        %     in2 = flip(in); % create a flipped in-grid (need row, column instead of column, row)
        paramsin = tif(in); % identify parameter values in the footprint
        if sum(sum(in)) > 0
            parameter_report(t) = nanmean(paramsin);
            range(t) = nanmax(paramsin) - nanmin(paramsin);
            SDev(t) = nanstd(paramsin);
        end
        clear xv yv in paramsin;
    end
elseif nargin == 2
    for k = 1:length(tif)
        %define the reference elevation data
        tif(k).x = tif(k).R.XWorldLimits(1)+0.5*tif(k).R.CellExtentInWorldX:tif(k).R.CellExtentInWorldX:tif(k).R.XWorldLimits(end)-0.5*tif(k).R.CellExtentInWorldX; % get a vector of x coords
        tif(k).y = tif(k).R.YWorldLimits(2)-0.5*tif(k).R.CellExtentInWorldY:-tif(k).R.CellExtentInWorldY:tif(k).R.YWorldLimits(1)+0.5*tif(k).R.CellExtentInWorldY; % get a vector of y coords
        [tif(k).xgrid, tif(k).ygrid] = meshgrid(tif(k).x, tif(k).y); % create grids of each of the x and y coords
        DEMdate(k) = tif(k).deciyear;
    end
    parameter_report = NaN([1, length(xs)]);
        
    %identify the reference elevation points in each ICESat2 footprint for
    %the dataset with the closest preceding date
    for t = 1:length(xs)
        xv = xs{t}; % bounding box x vector
        yv = ys{t}; % bounding box y vector
        
        %identify the appropriate terrain parameter dataset
        yr_ref = find(DEMdate<=T.date(t),1,'last');
        
        % data in the footprint
        in = inpolygon(tif(yr_ref).xgrid, tif(yr_ref).ygrid, xv, yv); % get logical array of in values
        paramsin = tif(yr_ref).z(in); % identify parameter values in the footprint
        if sum(sum(in)) > 0
            parameter_report(t) = nanmean(paramsin);
            range(t) = nanmax(paramsin) - nanmin(paramsin);
            SDev(t) = nanstd(paramsin);
        end
        clear xv yv in paramsin yr_ref;
    end
else
    disp('Non-compliant format for terrain parameter file');
end
params = parameter_report(:); % get the parameter column


