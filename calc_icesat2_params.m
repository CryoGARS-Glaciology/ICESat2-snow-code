function [params] = calc_icesat2_params(icesat2, tifs ,ATL0X)
% Function calc_icesat2_params extracts a specified terrain parameter from
% ICESat-2 footprints
% INPUTS: icesat2 = a csv file with icesat 2 elevations created using the
%                       h5 to csv jupyter notebook
%             tifs = the reference matrix or structure of the terrain
%                       parameter (slope, aspect, etc)
%           ATL0X = 6 or 8 to denote atl06 or alt08 respectively
% OUTPUTS: params = the terrain parameter from the tif averaged over the
%                   bounding box of the icesat2 footprint. Reported as a
%                   column vector matching the icesat2 input points

% last modified 03 March 2022 by Ellyn Enderlin (ellynenderlin@boisestate.edu)

[tif,Ref] = readgeoraster(tifs.name);
tif(tif==-9999) = NaN;

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
dates = datetime(T.time.Year,T.time.Month,T.time.Day);
[~,unique_refs] = unique(dates);
end_flag = zeros(size(norths,1),1);
end_flag(unique_refs) = 1; end_flag(unique_refs(unique_refs~=1)-1) = 1; end_flag(end) = 1;

%define the Reference elevation data
if isfield(Ref,'LatitudeLimits')
    [latgrid,longrid] = meshgrid(Ref.LongitudeLimits(1)+0.5*Ref.CellExtentInLongitude:Ref.CellExtentInLongitude:Ref.LongitudeLimits(2)-0.5*Ref.CellExtentInLongitude,...
        Ref.LatitudeLimits(2)-0.5*Ref.CellExtentInLatitude:-Ref.CellExtentInLatitude:Ref.LatitudeLimits(1)+0.5*Ref.CellExtentInLatitude);
    [xgrid, ygrid,~] = wgs2utm(latgrid,longrid);
else
    x = Ref.XWorldLimits(1)+0.5*Ref.CellExtentInWorldX:Ref.CellExtentInWorldX:Ref.XWorldLimits(2)-0.5*Ref.CellExtentInWorldX;
    if strcmp(Ref.ColumnsStartFrom,'north')
        y = Ref.YWorldLimits(2)-0.5*Ref.CellExtentInWorldY:-Ref.CellExtentInWorldY:Ref.YWorldLimits(1)+0.5*Ref.CellExtentInWorldY;
    else
        y = Ref.YWorldLimits(1)+0.5*Ref.CellExtentInWorldY:Ref.CellExtentInWorldY:Ref.YWorldLimits(2)-0.5*Ref.CellExtentInWorldY;
    end
    [xgrid, ygrid] = meshgrid(x, y); % create grids of each of the x and y coords
end

% calculates footprint corners
[xc,yc,theta] = ICESat2_FootprintCorners(norths,easts,ATL0X,end_flag);
parameter_report = NaN([1, size(xc,1)]);

%identify the reference elevation points in each ICESat2 footprint for
%the dataset with the closest preceding date

DEMdate = '20120826'; %only need to change this if the DTM is a geotiff

for k = 1:size(xc,1)
    xv = [xc(k,:) xc(k,1)]; % bounding box x vector
    yv = [yc(k,:) yc(k,1)]; % bounding box y vector

    %data in the footprint
    in = inpolygon(xgrid, ygrid, xv, yv); % get logical array of in values
  
    % data in the footprint
    paramsin = tif(in); % isolate parameter within the footprint
    if sum(sum(in)) > 0
        parameter_report(k) = nanmean(paramsin);
    else
        disp('Non-compliant format for terrain parameter file');
    end
    clear xv yv in paramsin yr_ref;
end
params = parameter_report(:); % get the parameter column


