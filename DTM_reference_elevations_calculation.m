%%% This code was writen by Karina Zikan and Ellyn Enderlin to calculate
%%% non-weighted mean, weighted mean, and weighted fitted refference
%%% elevations from a reference DTM and mean slope and aspect from slope 
%%% and aspect maps to compare to ICESat-2 ATL08 data
%%% (Use DTM_reference_elevations_calculation_sliderulecsv.m for ATL06 data from sliderule)
%%% << Adding mean vegitation is planned >>
%%% 
%%% SPECIFIED INPUTS:
%%%     DTM_path = path to the reference DTM on your computer
%%%     DTM_name = DTM file name
%%%     DTM_slope = Slope map file name
%%%     DTM_aspect = Aspect map file name
%%%     csv_path = path to the ICESat-2 datafiles on your computer
%%%     csv_name = name of ICESat-2 csv file
%%%     abbrev = site abriviation for file name
%%%     acronym = ICESat-2 product acronym
%%% OUTPUTS:
%%%     Reference_Elevations = csv datatable reporting the non-weighted
%%%         mean, std, weighted mean, and weighted fitted refference elevations,
%%%         mean slope, std slope, mean aspect, std aspect
%%%         
%%%
%%% Last updated: Jan 2023 by Karina Zikan


%% Inputs
clearvars; close all;
addpath(['./functions']) 

%DTM (be sure the path ends in a /)
DTM_path = 'DryCreek/DCEW_DEM/';
DTM_name = 'DryCreekBase1m_WGS84UTM11_WGS84.tif';
if contains(DTM_name,'.tif')
    DTM_date = '20120826'; %only need to change this if the DTM is a geotiff
end
% Slope
DTM_slope = 'DryCreekBase1m_WGS84UTM11_WGS84_slope.tif';
% Aspect
DTM_aspect = 'DryCreekBase1m_WGS84UTM11_WGS84_aspect.tif';


%csv (be sure the path ends in a /)
csv_path = '/Users/karinazikan/Documents/ICESat2-snow-code/DryCreek/';

%site abbreviation for file names
abbrev = 'DCEW';

%ICESat-2 product acronym
acronym = 'ATL08';
if acronym == 'ATL08'
    ATL0X = 8;
elseif acronym == 'ATL06'
    ATL0X = 6;
else
    error('acronym must be ATL06 or ATL08')
end 

%% Read in files
%days of year
modays_norm = [31 28 31 30 31 30 31 31 30 31 30 31];
cumdays_norm = cumsum(modays_norm); cumdays_norm = [0 cumdays_norm(1:11)];
modays_leap = [31 29 31 30 31 30 31 31 30 31 30 31];
cumdays_leap = cumsum(modays_leap); cumdays_leap = [0 cumdays_leap(1:11)];

%read in the snow-off reference elevation map
cd_to_DTM = ['cd ',DTM_path]; eval(cd_to_DTM);
if contains(DTM_name,'.tif')
    [DTM,Ref] = readgeoraster(DTM_name);
    DEMdate = DTM_date;
elseif contains(DTM_name,'.mat')
    load_DTM = ['load ',DTM_name]; eval(load_DTM);
    for i = 1:length(Z)
        DEMdate(i) = Z(i).deciyear;
    end
end
slope = readgeoraster(DTM_slope);
aspect = readgeoraster(DTM_aspect);

%identify the ICESat-2 data csv files
cd_to_csv = ['cd ',csv_path]; eval(cd_to_csv);
csvs = dir([acronym,'*.csv']); %if running more than once, rename original csvs to '*raw.csv' and change the search here to find files with that ending

%filter R2erence DTM elevations
% elevations(elevations < -10) = nan;
elevations = DTM;
elevations(elevations < -10) = nan; % throw out trash data
elevations(elevations > 10000) = nan; % more trash takeout

%load the ICESat-2 data
T = table; %create a table
icesat2 = [csv_path abbrev '-ICESat2-ATL08-params']; %compile the file name
file = readtable(icesat2); %read in files
T = [T; file];
%T = T(1:5,:); % ONLY FOR TESTING!!!!!!!!!!

% T = T([1:250],:);
zmod = T.Elevation(:); % save the median 'model' elevations (icesat-2 elevations)
zmodfit = T.Elevation_bestfit(:); % save the fitted 'model' elevations (icesat-2 elevations_bestfit)
zmodfit(isnan(zmod)) = NaN;
zstd = T.std; %save the standard deviation of the icesat-2 elevation estimates
easts = T.Easting(:); % pull out the easting values
norths = T.Northing(:); % pull out the northings
footwidth = 11; % approx. width of icesat2 shot footprint in meters

%identify the ends of each transect and flag them so that neighboring
%transects aren't used when constructing footprints (use beam variable & date)
dates = T.date;
[~,unique_refs] = unique([num2str(dates)],'rows');
end_flag = zeros(size(norths,1),1);
end_flag(unique_refs) = 1; end_flag(unique_refs(unique_refs~=1)-1) = 1; end_flag(end) = 1;


%% Calculating footprints for each data point
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

%% Calculate Reference Elevations, Slope, & Aspect
for r=1:length(zmod)  
    %identify the R2erence elevation points in each ICESat2 footprint
    xv = xc(r,[3:6 3]); % bounding box x vector
    yv = yc(r,[3:6 3]); % bounding box y vector

    %data in the footprint
    in = inpolygon(xgrid, ygrid, xv, yv); % get logical array of in values
    pointsinx = xgrid(in); % save x locations
    pointsiny = ygrid(in); % save y locations
    elevationsin = elevations(in); % save elevations
    slopesin = slope(in); % save slopes
    aspectsin = aspect(in); % save slopes
    
    %wieghted average
    dist = nan([1,length(pointsinx)])'; %initialize dist
    for a = 1:length(pointsinx)
        phi = atan2d((pointsiny(a)-norths(r)),(pointsinx(a)-easts(r)));
        dist(a)=abs(sqrt((pointsiny(a)-norths(r))^2+(pointsinx(a)-easts(r))^2)*sind(phi-theta(r))); %distance from the line in the center of the window  
    end
    maxdist = footwidth/2; % defining the maximum distance a point can be from the center icesat2 point
    w = 15/16*(1-(dist/maxdist).^2).^2; %bisqared kernel
    elevation_report_mean(r,:) = sum(w.*elevationsin)./sum(w); %weighted elevation estimate
    elevation_report_std(r,:) = std(elevationsin); %std of the elevations within the footprint

    %non wieghted average
    elevation_report_nw_mean(r,:) = nanmean(elevationsin); % non-wieghted elevations
    slope_mean(r,:) = nanmean(slopesin);
    slope_std(r,:) = std(slopesin);
    aspect_mean(r,:) = nanmean(aspectsin);
    aspect_std(r,:) = std(aspectsin);

    %weighted fit
    warning('off')
    p{1} = fit([pointsinx, pointsiny],elevationsin,'poly11','Weights',w); %fit linear polynomial
    p{2} = fit([pointsinx, pointsiny],elevationsin,'poly33','Weights',w); %fit cubic polynomial
    p{3} = fit([pointsinx, pointsiny],elevationsin,'poly44','Weights',w); %fit quadratic polynomial
    warning('on')
    for n=1:length(p) %loop through the degrees in d
        Em(n) = p{n}(easts(r),norths(r)); % Evaluate the fitted polynomial
        fitted = p{n}(pointsinx, pointsiny); % model elevation at each DEM location
        if n == 1 %calculatating corected linear midpoint elevation
            SlopeCorectedHight = fitted-elevationsin+Em(n);
            DistAlongWeight = 1/sqrt((pointsiny-norths(r)).^2+(pointsinx-easts(r)).^2-dist.^2);
            Em(n) = sum(SlopeCorectedHight.*DistAlongWeight,'all')/sum(DistAlongWeight,'all');
        end
        %RMSE(n) = sqrt(mean((Em(:,n)-elevationsin).^2)); % Calculate RMSE of fitted polynomial p
        fitstd(n) = std(fitted-elevationsin);
        fitmean(n) = mean(fitted-elevationsin);
    end
    a = find(fitstd==min(fitstd));
    if length(a) ~= 1
        fitmean = fitmean(a);
        a = find(fitmean==min(fitmean));
    end
    if length(a) ~= 1
        a = find(a==max(a));
    end
    order(r) = a;
    elevation_report_fitted(r,:) = Em(a);
end

%Write reference elevation table
E = table(elevation_report_nw_mean,elevation_report_mean,elevation_report_fitted,elevation_report_std,slope_mean,slope_std,aspect_mean,aspect_std);
writetable(E,[abbrev,'-ICESat2-',acronym,'-ref-elevations.csv']);


%% Sanity Checks
% % Surface fits
% figure;
% subplot(1,3,1);
% scatter3(pointsinx,pointsiny,elevationsin,'.','m')
% hold on
% plot(p{1})
% xlabel('Easting (km)')
% ylabel('Northing (km)')
% zlabel('Elevation (m)')
% colormap gray
% set(gca,'FontSize',16)
% label('Linear')
% 
% subplot(1,3,2);
% scatter3(pointsinx,pointsiny,elevationsin,'.','m')
% hold on
% plot(p{2})
% xlabel('Easting (km)')
% ylabel('Northing (km)')
% zlabel('Elevation (m)')
% colormap gray
% set(gca,'FontSize',16)
% label('Cubic')
% 
% subplot(1,3,3);
% scatter3(pointsinx,pointsiny,elevationsin,'.','m')
% hold on
% plot(p{3})
% xlabel('Easting (km)')
% ylabel('Northing (km)')
% zlabel('Elevation (m)')
% colormap gray
% set(gca,'FontSize',16)
% label('Quadratic')
% 
% % Distance and weighting check
% figure;
% plot3(pointsinx, pointsiny,dist,'.')
% xlabel('Easting (km)')
% ylabel('Northing (km)')
% zlabel('Distance from ICESat-2 track centerline (m)')
% 
% figure;
% plot3(pointsinx, pointsiny,w,'.')
% xlabel('Easting (km)')
% ylabel('Northing (km)')
% zlabel('Point weighting (m)')




