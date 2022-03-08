function [params] = calc_icesat2_params(icesat2, tif, R2)
% Function calc_icesat2_params extracts a specified terrain parameter from
% ICESat-2 footprints
% INPUTS: icesat2 = a csv file with icesat 2 elevations created using the
%                       h5 to csv jupyter notebook
%             tif = the reference matrix or structure of the terrain
%                       parameter (slope, aspect, etc)
%              R2 = the cell map reference
% OUTPUTS: params = the terrain parameter from the tif averaged over the
%                   bounding box of the icesat2 footprint. Reported as a
%                   column vector matching the icesat2 input points

% last modified 03 March 2022 by Ellyn Enderlin (ellynenderlin@boisestate.edu)

%specify ICESat-2 footprint width & length
footwidth = 11; % approx. width of icesat2 shot footprint in meters
if contains(icesat2(1,:), 'ATL08') % ATL08 commands
    default_length = 100; % approx. length of icesat2 shot footprint in meters
elseif contains(icesat2(1,:), 'ATL06') % ATL06 commands
    default_length = 40; % approx. length of icesat2 shot footprint in meters
end

%read the ICESat-2 data
T = readtable(icesat2);
easts = T.Easting(:); % pull out the easting values
norths = T.Northing(:); % pull out the northings

%identify the ends of each transect and flag them so that neighboring
%transects aren't used when constructing footprints (use beam variable & date)
beams = T.beam; dates = T.date;
[~,unique_refs] = unique([num2str(dates),string(beams)],'rows');
end_flag = zeros(size(norths,1),1);
end_flag(unique_refs) = 1; end_flag(unique_refs(unique_refs~=1)-1) = 1; end_flag(end) = 1;

% initialize matrix for RGT orientations
theta = NaN(size(norths,1),2);

%create polygons of ICESat-2 footprints
for r = 1:size(theta,1)
    %calculate the footprint orientation
    if r == 1  || end_flag(r)-end_flag(r-1) == 0 %only angle pointed forwards
        theta(r,1) = atan2d((norths(r+1)-norths(r)),(easts(r+1)-easts(r))); theta(r,2) = theta(r,1);
        footlength = default_length;
    elseif r == size(theta,1) || end_flag(r)-end_flag(r+1) == 0 %only angle pointed backwards
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

%extract data from the input geotiff if over static terrain (ATL08) or from
%the geotiff timeseries if over a temporally-evolving surface (ATL06)
if nargin == 3
    %define the reference elevation data
    x = R2.XWorldLimits(1)+0.5*R2.CellExtentInWorldX:R2.CellExtentInWorldX:R2.XWorldLimits(end)-0.5*R2.CellExtentInWorldX; % get a vector of x coords
    y = R2.YWorldLimits(2)-0.5*R2.CellExtentInWorldY:-R2.CellExtentInWorldY:R2.YWorldLimits(1)+0.5*R2.CellExtentInWorldY; % get a vector of y coords
    [xgrid, ygrid] = meshgrid(x, y); % create grids of each of the x and y coords
    parameter_report = NaN([1,length(xc)]);
    
    %identify the reference elevation points in each ICESat2 footprint
    for k = 1:length(xc)
        xv = [xc(k,:) xc(k,1)]; % bounding box x vector
        yv = [yc(k,:) yc(k,1)]; % bounding box y vector
        
        %data in the footprint
        in = inpolygon(xgrid, ygrid, xv, yv); % get logical array of in values
        paramsin = tif(in); % isolate parameter within the footprint
        if sum(sum(in)) > 0
            parameter_report(k) = nanmean(paramsin);
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
    parameter_report = NaN([1, size(xc,1)]);
    
    %identify the reference elevation points in each ICESat2 footprint for
    %the dataset with the closest preceding date
    for k = 1:size(xc,1)
        xv = [xc(k,:) xc(k,1)]; % bounding box x vector
        yv = [yc(k,:) yc(k,1)]; % bounding box y vector
        
        %identify the appropriate terrain parameter dataset
        yr_ref = find(DEMdate<=T.date(k),1,'last');
        
        % data in the footprint
        in = inpolygon(tif(yr_ref).xgrid, tif(yr_ref).ygrid, xv, yv); % get logical array of in values
        paramsin = tif(yr_ref).z(in); % isolate parameter within the footprint
        if sum(sum(in)) > 0
            parameter_report(k) = nanmean(paramsin);
        end
        clear xv yv in paramsin yr_ref;
    end
else
    disp('Non-compliant format for terrain parameter file');
end
params = parameter_report(:); % get the parameter column


