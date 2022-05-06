function [params] = calc_icesat2_params(icesat2, tif, R2)
% Function CALC_ICESAT2_PARAMS extracts a specified terrain parameter from
% ICESat-2 segments.
% INPUTS: icesat2 = csv file(s) with ICESat-2 data saved as column vectors
%      elevations = the reference elevation matrix 
%              R2 = the cell map reference for the reference DTM
% OUTPUTS: params = the terrain parameter from the tif averaged over the
%                   bounding box of the icesat2 footprint. Reported as a
%                   column vector matching the icesat2 input points

% last modified 03 March 2022 by Ellyn Enderlin (ellynenderlin@boisestate.edu)

%specify ICESat-2 footprint width & length
footwidth = 11; % approx. width of icesat2 shot footprint in meters
if contains(icesat2(1,:), 'ATL08') % ATL08 commands
    default_length = 100; % approx. length of icesat2 shot footprint in meters
    ATL0X = 8; %dataset flag for footprint delineation
elseif contains(icesat2(1,:), 'ATL06') % ATL06 commands
    default_length = 40; % approx. length of icesat2 shot footprint in meters
    ATL0X = 6; %dataset flag for footprint delineation
end

%read the ICESat-2 data
T = readtable(icesat2);
easts = T.Easting(:); % pull out the easting values
norths = T.Northing(:); % pull out the northings

%identify the ends of each transect and flag them so that neighboring
%transects aren't used when constructing segments (use beam variable & date)
beams = T.beam; dates = T.date;
[~,unique_refs] = unique([num2str(dates),string(beams)],'rows');
end_flag = zeros(size(norths,1),1);
end_flag(unique_refs) = 1; end_flag(unique_refs(unique_refs~=1)-1) = 1; end_flag(end) = 1;

%create polygons of ICESat-2 segments
[xc,yc,~] = ICESat2_FootprintCorners(norths,easts,ATL0X,end_flag);

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


