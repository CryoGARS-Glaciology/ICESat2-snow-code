%%% This code was writen by Karina Zikan to calculate snow depth statistics
%%% from ICESat-2 and make visualization plots
%%%
%%% SPECIFIED INPUTS:
%%%     icesat2 = path to ICESat-2 file (output from batch_icesat2_terrain_parameter_comparison.m)
%%%     ref_elevations = path to reference elevations file (output from DTM_reference_elevations_calculation.m)
%%%     abbrev = site abriviation for file name
%%%     snowcover = set to snow-on ('snowon') or snow-off ('snowoff') conditions
%%% OUTPUTS:
%%%
%%%
%%% Last updated: November 2022 by Karina Zikan

%% Inputs
clearvars;
addpath(['./functions'])
addpath(['/Users/karinazikan/Documents/cmocean'])

%Folder path
folderpath = '/Users/karinazikan/Documents/ICESat2-snow-code/RCEW/';
%site abbreviation for file names
abbrev = 'RCEW';

% DEM path
DTM_name = [folderpath 'RCEW_DEM/RCEW_1m_WGS84UTM11_WGS84.tif'];

%%
%File paths
icesat2_atl08 = [folderpath abbrev '-ICESat2-ATL08-params'];
ref_elevations_atl08 = [folderpath abbrev '-ICESat2-ATL08-ref-elevations'];

icesat2_atl06 = [folderpath abbrev '-ICESat2-ATL06sr-params'];
ref_elevations_atl06 = [folderpath abbrev '-ICESat2-ATL06sr-ref-elevations'];

icesat2_atl06_class = [folderpath abbrev '-ICESat2-ATL06sr-atl08class-params'];
ref_elevations_atl06_class = [folderpath abbrev '-ICESat2-ATL06sr-atl08class-ref-elevations'];

%set colors
colors{1} = cmocean('-dense',6);
colors{2} = cmocean('-algae',5);
colors{3} = cmocean('ice',5);

%% Load data
%load the reference elevation data
E_08 = readtable(ref_elevations_atl08);
E_06 = readtable(ref_elevations_atl06);
E_06class = readtable(ref_elevations_atl06_class);

%load the ICESat-2 data
I_08 = readtable(icesat2_atl08); %read in files
I_06 = readtable(icesat2_atl06);
I_06class = readtable(icesat2_atl06_class);

footwidth = 11; % approx. width of icesat2 shot footprint in meters

%DEM 
[DTM,Ref] = readgeoraster(DTM_name);
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
DTM(DTM>(2.5*10^3)) = NaN;
x = x./(10^3); % to km
y = y./(10^3);% to km
%y = flip(y);
DTM(DTM<0) = NaN;

%% Make snow-on and snow-off arrays ({1} = ATL08, {2} = ATL06, {3} = ATL06 w/ ATL08 classification)

%ATL08
% season (1=winter(Jan,Feb,Mar),2=spring(Apr,May,Jun),3=summer(Jul,Aug,Sep),4=fall(Oct,Nov,Dec))
ix_off = find(I_08.season == 3);
ix_on = find(I_08.season == 1);
I_off{1} = I_08(ix_off,:);
E_off{1} = E_08(ix_off,:);
I_on{1} = I_08(ix_on,:);
E_on{1} = E_08(ix_on,:);

%ATl06
ix_off = find(I_06.time.Month >=6 & I_06.time.Month <= 9);
ix_on = find(I_06.time.Month <=2 | I_06.time.Month >=12);
I_off{2} = I_06(ix_off,:);
E_off{2} = E_06(ix_off,:);
I_on{2} = I_06(ix_on,:);
E_on{2} = E_06(ix_on,:);

%ATl06 atl08class
ix_off = find(I_06class.time.Month >=6 & I_06class.time.Month <= 9);
ix_on = find(I_06class.time.Month <=2 | I_06class.time.Month >=12);
I_off{3} = I_06class(ix_off,:);
E_off{3} = E_06class(ix_off,:);
I_on{3} = I_06class(ix_on,:);
E_on{3} = E_06class(ix_on,:);

%% Calculate Snow-on and Snow-off residuals & snow depth (snow-on residuals shifted by meadian Snow-off residual)
for i = 1:3
    if i == 1
%         Residuals_off{i} =  I_off{i}.Elevation_bestfit - E_off{i}.elevation_report_fitted;
%        Residuals_on{i} =   I_on{i}.Elevation_bestfit - E_on{i}.elevation_report_fitted;
        Residuals_off{i} =  I_off{i}.Elevation - E_off{i}.elevation_report_fitted;
        Residuals_on{i} =   I_on{i}.Elevation - E_on{i}.elevation_report_fitted;
    else
        Residuals_off{i} =  I_off{i}.h_mean - E_off{i}.elevation_report_nw_mean;
        Residuals_on{i} =  I_on{i}.h_mean - E_on{i}.elevation_report_nw_mean;
    end
    Residuals_off{i}(Residuals_off{i} > 80) = NaN; Residuals_off{i}(Residuals_off{i} < -80) = NaN; %remove extreme outliers
    Residuals_on{i}(Residuals_on{i} > 80) = NaN; Residuals_on{i}(Residuals_on{i} < -80) = NaN; %remove extreme outliers
    SnowDepth{i} = Residuals_on{i} - median(Residuals_off{i},'omitnan');
end

%% Figures
fig1 = figure(1); clf % histgrams with non-parametric pdfs
binwidth = 0.2;
subplot(3,1,1); hold on %snow off residuals
for i = 1:3
    h(i) = histogram(Residuals_off{i},'Normalization','pdf');  h(i).FaceAlpha = .15; h(i).BinWidth = binwidth; h(i).FaceColor = colors{i}(3,:);  h(i).EdgeColor = 'k';
    pd = fitdist(Residuals_off{i},'kernel','Kernel','normal');
    fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2,'Color',colors{i}(3,:));
end
%plot([0,0],[0,.8], 'linewidth', 2, 'Color','k') % plot reference 0 line
set(gca,'fontsize',16,'xlim',[-4 2]);
set(gcf,'position',[50 50 800 400]);
xlabel('Snow off Vertical offset (m)'); ylabel('Probability density');
txt = {['N-08 = ' num2str(length(SnowDepth{1})-sum(isnan(SnowDepth{1})))],['N-06 = ' num2str(length(SnowDepth{2})-sum(isnan(SnowDepth{2})))],['N-06class = ' num2str(length(SnowDepth{3})-sum(isnan(SnowDepth{3})))]};
text(-3,.5,txt);

subplot(3,1,2); hold on %snow on residuals
for i = 1:3
    h(i) = histogram(Residuals_on{i},'Normalization','pdf');  h(i).FaceAlpha = .15; h(i).BinWidth = binwidth; h(i).FaceColor = colors{i}(3,:);  h(i).EdgeColor = 'k';
    pd = fitdist(Residuals_on{i},'kernel','Kernel','normal');
    fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2,'Color',colors{i}(3,:));
end
%plot([0,0],[0,.8], 'linewidth', 2, 'Color','k') % plot reference 0 line
set(gca,'fontsize',16,'xlim',[-4 2]);
set(gcf,'position',[50 50 800 400]);
xlabel('Snow On Vertical offset (m)'); ylabel('Probability density');
legend('ATL08','ATL08 pdf','ATL06sr','ATL06sr pdf','ATL06sr ATL08 class','ATL06sr ATL08 class pdf','Location','best');

subplot(3,1,3); hold on %snow depth
for i = 1:3
    h(i) = histogram(SnowDepth{i},'Normalization','pdf');  h(i).FaceAlpha = .15; h(i).BinWidth = binwidth; h(i).FaceColor = colors{i}(3,:);  h(i).EdgeColor = 'k';
    pd = fitdist(SnowDepth{i},'kernel','Kernel','normal');
    fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2,'Color',colors{i}(3,:));
end
%plot([0,0],[0,.8], 'linewidth', 2, 'Color','k') % plot reference 0 line
set(gca,'fontsize',16,'xlim',[-4 2]);
set(gcf,'position',[50 50 800 400]);
xlabel('Snow Depth(m)'); ylabel('Probability density');
txt = {['Median-08 = ' num2str(median(SnowDepth{1},'omitnan'))],
    ['Median-06 = ' num2str(median(SnowDepth{2},'omitnan'))],
    ['Median-06class = ' num2str(median(SnowDepth{3},'omitnan'))]};
text(-3,.5,txt);

%-------------------------------------------------------------------------

fig2 = figure(2); clf %non-parametric pdfs
binwidth = 0.2;
subplot(3,1,1); hold on %snow off residuals
for i = 1:3
    pd = fitdist(Residuals_off{i},'kernel','Kernel','normal');
    fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2,'Color',colors{i}(3,:));
end
%plot([0,0],[0,.8], 'linewidth', 2, 'Color','k') % plot reference 0 line
set(gca,'fontsize',16,'xlim',[-3 2]);
set(gcf,'position',[50 50 800 400]);
xlabel('Snow off Vertical offset (m)'); ylabel('Probability density');
txt = {['N-08 = ' num2str(length(Residuals_off{1})-sum(isnan(Residuals_off{1})))],['N-06 = ' num2str(length(Residuals_off{2})-sum(isnan(Residuals_off{2})))],['N-06class = ' num2str(length(Residuals_off{3})-sum(isnan(Residuals_off{3})))]};
text(-2.9,.8,txt);

subplot(3,1,2); hold on %snow on residuals
for i = 1:3
    pd = fitdist(Residuals_on{i},'kernel','Kernel','normal');
    fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2,'Color',colors{i}(3,:));
end
%plot([0,0],[0,.8], 'linewidth', 2, 'Color','k') % plot reference 0 line
set(gca,'fontsize',16,'xlim',[-3 2]);
set(gcf,'position',[50 50 800 400]);
xlabel('Snow On Vertical offset (m)'); ylabel('Probability density');
legend('ATL08','ATL06sr','ATL06sr w/ 08classes','Location','best');
txt = {['N-08 = ' num2str(length(Residuals_on{1})-sum(isnan(Residuals_on{1})))],['N-06 = ' num2str(length(Residuals_on{2})-sum(isnan(Residuals_on{2})))],['N-06class = ' num2str(length(Residuals_on{3})-sum(isnan(Residuals_on{3})))]};
text(-2.9,.8,txt);

subplot(3,1,3); hold on %snow depth
for i = 1:3
    pd = fitdist(SnowDepth{i},'kernel','Kernel','normal');
    fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2,'Color',colors{i}(3,:));
end
%plot([0,0],[0,.8], 'linewidth', 2, 'Color','k') % plot reference 0 line
set(gca,'fontsize',16,'xlim',[-3 2]);
set(gcf,'position',[50 50 800 400]);
xlabel('Snow Depth(m)'); ylabel('Probability density');
txt = {['Median-08 = ' num2str(median(SnowDepth{1},'omitnan'))],
    ['Median-06 = ' num2str(median(SnowDepth{2},'omitnan'))],
    ['Median-06class = ' num2str(median(SnowDepth{3},'omitnan'))]};
text(-2.9,.8,txt);

fig3 = figure(3); clf; hold on
for i = 1:3
    subplot(3,1,i); hold on
    pd = fitdist([Residuals_off{i}-median(Residuals_off{i}, 'omitnan')],'kernel','Kernel','normal');
    fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2,'Color',colors{i}(2,:));
    pd = fitdist([Residuals_on{i}-median(Residuals_off{i}, 'omitnan')],'kernel','Kernel','normal');
    fplot(@(x) pdf(pd,x),[-10 8], 'Linewidth', 2,'Color',colors{i}(4,:));
    set(gca,'fontsize',16,'xlim',[-3 3]);
    set(gcf,'position',[50 50 800 400]);
    xlabel('Vertical offset (m)'); ylabel('Probability density');
    legend('Summer','Winter');
end
subplot(3,1,1)
title('ATL08');
subplot(3,1,2)
title('ATL06');
subplot(3,1,3)
title('ATL06 with ATL08 photon classification');
hold off

%-------------------------------------------------------------------------
%% Grouped data plots
clear SnowDepthTable
% ATL06 classified
% group elevation
Residuals =  I_06class.h_mean - E_06class.elevation_report_nw_mean;
Residuals(Residuals > 80) = NaN; Residuals(Residuals < -80) = NaN; %remove extreme outliers
SnowDepthAll = Residuals - median(Residuals_off{3},'omitnan');

E_06class.elevation_report_nw_mean(isnan(SnowDepthAll)) = NaN;

figure(5);
h = histogram(E_06class.elevation_report_nw_mean,30);
elev_binwidth = h.BinWidth; elev_binedges = h.BinEdges;
clear h;
for i = 2:length(elev_binedges)
    bins{i-1} = num2str(elev_binedges(i));
end
groupElev = discretize(E_06class.elevation_report_nw_mean, elev_binedges,'categorical',bins);

% Make table that can be grouped by various parameters
SnowDepthTable(:,1) = table([datetime(I_06class.time.Year,I_06class.time.Month,I_06class.time.Day)]); %Date
SnowDepthTable(:,2) = table(groupElev);  %binned elevations
SnowDepthTable(:,3) = table(SnowDepthAll); %ATL06_class Snow Depth
SnowDepthTable = renamevars(SnowDepthTable,["Var1","Var2","Var3"],["Date","Elevation","ATL06_class Snow Depth"]);

% Average snow depth by elevation by date
SnowDepthMedTable = varfun(@(x)median(x,'omitnan'),SnowDepthTable,'GroupingVariables',{'Date','Elevation'});
SnowDepthMedDatesTable = varfun(@(x)median(x,'omitnan'),SnowDepthTable,'InputVariables','ATL06_class Snow Depth','GroupingVariables','Date');
stdT = varfun(@std,SnowDepthTable,'GroupingVariables',{'Date','Elevation'});
SnowDepthMedTable(:,5) = stdT(:,4); 
SnowDepthMedTable = renamevars(SnowDepthMedTable,["Fun_ATL06_class Snow Depth","Var5"],["median","std"]);
dates = unique(SnowDepthMedTable.Date);

% Make grouped snow depths array
for k = 1:length(dates)
    for i = 2:length(elev_binedges)
        ix = SnowDepthMedTable.Date == dates(k); %index by date
        T = SnowDepthMedTable(ix,:);
        if sum(ismember(num2str(elev_binedges(i)),T.Elevation)) ~= 0
            ix = find(T.Elevation == num2str(elev_binedges(i))); %index by elevation
            SnowDepthMedian = T.median(ix,:);
        else
            SnowDepthMedian = NaN;
        end
        SnowDepthMedArray(k,[i-1]) = SnowDepthMedian;
    end
end

%% Plot grouped snow depths
fig4 = figure(4); clf
datenum = convertTo(dates,'yyyymmdd');
imagesc(elev_binedges([2:length(elev_binedges)]),datenum,SnowDepthMedArray,[-3 3])
cmap = cmocean('-balance'); cmap = [ 0 0 0 ; cmap ];
colormap(cmap); c = colorbar;
set(gca,'fontsize',16);
xlabel('Elevation')
a=(max(datenum)-min(datenum))/length(datenum);
ticks = min(datenum):a:max(datenum);
yticks(ticks)
yticklabels(string(dates,'MM-yyyy'))
%-------------------------------------------------------------------------
% %% Single Track plots
% % ATL06 w/ atl08 classification
% date = '26-Jan-2020'; %pick track date
% dates = datetime(I_on{3}.time.Year,I_on{3}.time.Month,I_on{3}.time.Day);
% ix = find(dates == date);
% I_on_track = I_on{3}(ix,:);
% E_on_track = E_on{3}(ix,:);
% SnowDepth_track = SnowDepth{3}(ix,:);
% 
% Nstart=1;
% Nend=40;
% 
% fig5 = figure(5); clf
% subplot(2,1,1); hold on
% scatter([I_on_track.Easting(Nstart:Nend)./(10^3)],E_on_track.elevation_report_nw_mean(Nstart:Nend),'.')
% scatter([I_on_track.Easting(Nstart:Nend)./(10^3)],[I_on_track.h_mean(Nstart:Nend)-median(Residuals_off{3}, 'omitnan')],'.')
% xlabel('Easting (km)'); ylabel('Elevation (m)');
% legend('Snow Free DEM','ICESat-2')
% set(gca,'fontsize',16);
% set(gcf,'position',[50 50 800 400]);
% 
% subplot(2,1,2); hold on
% stem([I_on_track.Easting(Nstart:Nend)./(10^3)], SnowDepth_track(Nstart:Nend),'Linewidth', 2)
% xlabel('Easting (km)'); ylabel('Snow Depth (m)');
% %legend('Snow Depth')
% set(gca,'fontsize',16);
% set(gcf,'position',[50 50 800 400]);
% 
% %%
% % plot map + track
% fig5 = figure(5);
% imagesc(x,y,DTM) 
% daspect([1 1 1])
% colormap(cmocean('grey'))
% hold on
% scatter([I_on_track.Easting./(10^3)],[I_on_track.Northing./(10^3)],[],'g','.')
% xlabel('Easting (km)')
% ylabel('Northing (km)')
% set(gca,'fontsize',16);
% 
% 
