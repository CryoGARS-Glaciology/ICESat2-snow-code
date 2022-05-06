%checking ICESat2 footprint mapping onto DEMs
clear all; close all;
addpath(['./functions']) 

ATL0X = 8;

boxW = 11;
boxL = 100;

csv_path = '/Users/karinazikan/Documents/ICESat2-snow-code/';

%load the ICESat-2 data
T = table; %create a table
% for i = 1:length(csvs)
%     icesat2 = [csv_path,'RCEW-ICESat2-ATL08-params']; %compile the file name
%     file = readtable(icesat2); %read in files
%     T = [T; file]; %combine tables
% end
icesat2 = [csv_path,'RCEW-ICESat2-ATL08-params']; %compile the file name
file = readtable(icesat2); %read in files
T = [T; file];

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

% % for ATL08 files only use snow-free data (brightness flag == 0)
% if contains(icesat2, 'ATL08') % ATL08 commands
%     bright = T.Brightness_Flag;
%     ib = find(bright == 0);
%     easts = easts(ib);
%     norths = norths(ib);
%     zmod = zmod(ib);
%     zmodfit = zmodfit(ib);
% end

cd '/Users/karinazikan/Documents/ICESat2-snow-code'


% calculates footprint corners
[xc,yc,theta] = ICESat2_FootprintCorners(norths,easts,ATL0X,end_flag);
xc = xc(unique_refs,:);
yc = yc(unique_refs,:);
theta = theta(unique_refs,:);

%loop through each unidque date to check footprints to make sure that the
%code applies universally
for j = 1:length(unique_refs)
    xv = xc(j,[3:6 3]); % bounding box x vector
    yv = yc(j,[3:6 3]); % bounding box y vector


    %plot
    figure; set(gcf,'position',[500 50 400 1200]);
    pl(1) = plot(T.Easting(unique_refs(j):unique_refs(j)+2),T.Northing(unique_refs(j):unique_refs(j)+2),'-xk','linewidth',1.5); axis xy equal; hold on;
    plot(xc(j,3),yc(j,3),'sm','markersize',12); plot(xc(j,4),yc(j,4),'dm','markersize',12);
    plot(xc(j,5),yc(j,5),'sc','markersize',12); plot(xc(j,6),yc(j,6),'dc','markersize',12);
    pl(2) = plot(xv,...
        yv,'--b','linewidth',2); hold on;
    
    %OLD CALCULATIONS
    theta = abs(atan((T.Northing(unique_refs(j)+2) - T.Northing(unique_refs(j)+1))/(T.Easting(unique_refs(j)+2) - T.Easting(unique_refs(j)+1))));
    % get the x and y vectors to form the polygon
    xpoly(1) = T.Easting(unique_refs(j)) + (boxW/2) - (boxW/2)*cos((pi/2) - theta); % calculate the 4 corners in the x direction
    xpoly(2) = T.Easting(unique_refs(j)) + (boxW/2) + (boxW/2)*cos((pi/2) - theta);
    xpoly(3) = T.Easting(unique_refs(j)) - (boxW/2) + (boxW/2)*cos((pi/2) - theta);
    xpoly(4) = T.Easting(unique_refs(j)) - (boxW/2) - (boxW/2)*cos((pi/2) - theta);
    ypoly(1) = T.Northing(unique_refs(j)) - (boxL/2) - (boxW/2)*cos((pi/2) - theta); % calculate the 4 corners in the y direction
    ypoly(2) = T.Northing(unique_refs(j)) - (boxL/2) + (boxW/2)*cos((pi/2) - theta);
    ypoly(3) = T.Northing(unique_refs(j)) + (boxL/2) + (boxW/2)*cos((pi/2) - theta);
    ypoly(4) = T.Northing(unique_refs(j)) + (boxL/2) - (boxW/2)*cos((pi/2) - theta);
    pl(3) = plot([xpoly, xpoly(1)],[ypoly, ypoly(1)],...
        ':','color',[1 0.5 0],'linewidth',2); 
    set(gca,'xlim',[min([xc(j) xpoly]) max([xc(j) xpoly])],'ylim',[min([yc(j) ypoly]) max([yc(j) ypoly])]);
    
    %add a legend to the plot
    leg = legend(pl,'RGT','new box','old box'); set(leg,'location','eastoutside');
    title([num2str(j),') ',num2str(unique_refs(j))]);
    
    %manually confirm accuracy of the box
    answer = questdlg('Does the new box look correct?',...
        'Box Check',...
        'Yes','No','No'); %third option is the default
    switch answer
        case 'Yes'
            box_flag = 1;
        case 'No'
            box_flag = 0;
    end
    disp('moving on...'); close all;
end