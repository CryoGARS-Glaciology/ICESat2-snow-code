function rmsez = coregister_icesat2(icesat2, elevations, R2, A)
% Function COREGISTER_ICESAT2 calculates the root mean square error in
% elevation for all the ICESat-2 segments in a track, using the mean
% elevation from each segment (T.Elevation) and the mean of the reference
% elevation dataset for the same footprint.
% INPUTS: icesat2 = csv file(s) with ICESat-2 data saved as column vectors
%      elevations = the reference elevation matrix 
%              R2 = the cell map reference for the reference DTM
%               A = a [2 1] vector that serves as the spatial offsets in
%                       the x and y directions (meters)
% OUTPUTS:  rmsez = the root mean squared difference between the icesat-2
%                       elevations and their corresponding (offset) DTM 
%                       elevations

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

%filter reference DTM elevations
% elevations(elevations < -10) = nan;
elevations(elevations < -10) = nan; % throw out trash data
elevations(elevations > 10000) = nan; % more trash takeout

%load the ICESat-2 data
T = []; end_flag = [];
for j = 1:size(icesat2,1)
    t = readtable(icesat2(j,:)); T = [T; t];
    end_flag = [end_flag; 1; zeros(length(t.Elevation)-2,1); 1];
    clear t;
end
zmod = T.Elevation(:); % save the 'model' elevations (icesat-2 elevations)
easts = T.Easting(:); % pull out the easting values
norths = T.Northing(:); % pull out the northings

% for ATL08 files only use snow-free data (brightness flag == 0)
if contains(icesat2(1,:), 'ATL08') % ATL08 commands
    bright = T.Brightness_Flag;
    ib = find(bright == 0);
    easts = easts(ib);
    norths = norths(ib);
    zmod = zmod(ib);
    end_flag = end_flag(ib);
end
end_flag(end) = 1;

%create polygons of ICESat-2 segments
[xc,yc,~] = ICESat2_FootprintCorners(norths,easts,ATL0X,end_flag);

% %plot (uncomment if you want to quality check
% if A(1) == 0 && A(2) == 0
%     figure; set(gcf,'position',[500 50 400 1200]);
%     pl(1) = plot(easts,norths,'.k','linewidth',1.5); axis xy equal; hold on;
%     for r = 1:length(theta)
%         pl(2) = plot([xc(r,:), xc(r,1)],[yc(r,:), yc(r,1)],'--b','linewidth',2); hold on;
%     end
%     leg = legend(pl,'RGT','segments'); set(leg,'location','eastoutside');
%     answer = questdlg('Do the segments look correct?',...
%         'Box Check',...
%         'Yes','No','No'); %third option is the default
%     switch answer
%         case 'Yes'
%             disp('moving on with coregistration...');
%         case 'No'
%             disp('FIX FOOTPRINT TRIG!'); return
%     end
% end

%define the reference elevation data
x = R2.XWorldLimits(1)+0.5*R2.CellExtentInWorldX:R2.CellExtentInWorldX:R2.XWorldLimits(2)-0.5*R2.CellExtentInWorldX; 
if strcmp(R2.ColumnsStartFrom,'north')
    y = R2.YWorldLimits(2)-0.5*R2.CellExtentInWorldY:-R2.CellExtentInWorldY:R2.YWorldLimits(1)+0.5*R2.CellExtentInWorldY; 
else
    y = R2.YWorldLimits(1)+0.5*R2.CellExtentInWorldY:R2.CellExtentInWorldY:R2.YWorldLimits(2)-0.5*R2.CellExtentInWorldY; 
end
[xgrid, ygrid] = meshgrid(x, y); % create grids of each of the x and y coords
elevation_report = zeros([1, length(xc)]);

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
