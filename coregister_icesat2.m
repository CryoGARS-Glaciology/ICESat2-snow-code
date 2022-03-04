function rmsez = coregister_icesat2(icesat2, elevations, R2, A)
% Function COREGISTER_ICESAT2 coregisters icesat-2 data with a corresponding digital
% terrain model 
% INPUTS: icesat2 = csv file(s) with icesat 2 elevations created using the
%                       h5 to csv jupyter notebook
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
elseif contains(icesat2(1,:), 'ATL06') % ATL06 commands
    default_length = 40; % approx. length of icesat2 shot footprint in meters
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
end

% % initialize empty matrices
theta = NaN(length(norths),2);

%create polygons of ICESat-2 footprints
for r = 1:length(theta)
    %calculate the footprint orientation
    if r == 1  || end_flag(r)-end_flag(r-1) == 0 %only angle pointed forwards
        theta(r,1) = atan2d((norths(r+1)-norths(r)),(easts(r+1)-easts(r))); theta(r,2) = theta(r,1);
        footlength = default_length;
    elseif r == length(theta) || end_flag(r)-end_flag(r+1) == 0 %only angle pointed backwards
        theta(r,1) = atan2d((norths(r)-norths(r-1)),(easts(r)-easts(r-1))); theta(r,2) = theta(r,1);
        footlength = default_length;
    else %calculate angles in each direction
        theta(r,1) = atan2d((norths(r)-norths(r-1)),(easts(r)-easts(r-1)));
        theta(r,2) = atan2d((norths(r+1)-norths(r)),(easts(r+1)-easts(r)));
        footlength = sqrt((norths(r+1)-norths(r-1)).^2 + (easts(r+1)-easts(r-1)).^2)/2;
    end
    
    %find box edges along the RGT
    back_x = A(1)+easts(r)-(footlength/2)*cosd(theta(r,1)); back_y = A(2)+norths(r)-(footlength/2)*sind(theta(r,1));
    front_x = A(1)+easts(r)+(footlength/2)*cosd(theta(r,2)); front_y = A(2)+norths(r)+(footlength/2)*sind(theta(r,2));
    
    %find box edges perpendicular to the centroid
    xc(r,1) = A(1)+easts(r)+(footwidth/2)*cosd(nanmean(theta(r,:))+90); yc(r,1) = A(2)+norths(r)+(footwidth/2)*sind(nanmean(theta(r,:))+90);
    xc(r,4) = A(1)+easts(r)+(footwidth/2)*cosd(nanmean(theta(r,:))-90); yc(r,4) = A(2)+norths(r)+(footwidth/2)*sind(nanmean(theta(r,:))-90);
    
    %solve for corner coordinates
    xc(r,2) = back_x+(footwidth/2)*cosd(theta(r,1)+90); yc(r,2) = back_y+(footwidth/2)*sind(theta(r,1)+90);
    xc(r,3) = back_x+(footwidth/2)*cosd(theta(r,1)-90); yc(r,3) = back_y+(footwidth/2)*sind(theta(r,1)-90);
    xc(r,5) = front_x+(footwidth/2)*cosd(theta(r,2)-90); yc(r,5) = front_y+(footwidth/2)*sind(theta(r,2)-90);
    xc(r,6) = front_x+(footwidth/2)*cosd(theta(r,2)+90); yc(r,6) = front_y+(footwidth/2)*sind(theta(r,2)+90);
    clear back_* front_*;
end

% %plot (uncomment if you want to quality check
% if A(1) == 0 && A(2) == 0
%     figure; set(gcf,'position',[500 50 400 1200]);
%     pl(1) = plot(easts,norths,'.k','linewidth',1.5); axis xy equal; hold on;
%     for r = 1:length(theta)
%         pl(2) = plot([xc(r,:), xc(r,1)],[yc(r,:), yc(r,1)],'--b','linewidth',2); hold on;
%     end
%     leg = legend(pl,'RGT','footprints'); set(leg,'location','eastoutside');
%     answer = questdlg('Do the footprints look correct?',...
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
