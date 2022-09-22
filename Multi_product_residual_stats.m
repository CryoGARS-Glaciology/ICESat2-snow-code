%%% This code was writen by Karina Zikan to statistically compare elevations
%%% from multiple ICESat-2 products to reference elevations calculated by
%%% DTM_reference_elevations_calculation.m
%%%
%%% SPECIFIED INPUTS:
%%%     icesat2 = path to ICESat-2 file (output from batch_icesat2_terrain_parameter_comparison.m)
%%%     ref_elevations = path to reference elevations file (output from DTM_reference_elevations_calculation.m)
%%%     abbrev = site abriviation for file name
%%%     snowcover = set to snow-on ('snowon') or snow-off ('snowoff') conditions
%%% OUTPUTS:
%%%
%%%
%%% Last updated: September 2022 by Karina Zikan

%% Inputs
clearvars;
addpath(['./functions'])
addpath(['/Users/karinazikan/Documents/cmocean'])

%File paths
% icesat2_atl08 = '/Users/karinazikan/Documents/ICESat2-snow-code/RCEW/RCEW-ICESat2-ATL08-params';
% ref_elevations_atl08 = '/Users/karinazikan/Documents/ICESat2-snow-code/RCEW/RCEW-ICESat2-ATL08-ref-elevations';

icesat2_atl06 = '/Users/karinazikan/Documents/ICESat2-snow-code/DryCreek/DCEW-ICESat2-ATL06sr-params';
ref_elevations_atl06 = '/Users/karinazikan/Documents/ICESat2-snow-code/DryCreek/DCEW-ICESat2-ATL06sr-ref-elevations';

icesat2_atl06_class = '/Users/karinazikan/Documents/ICESat2-snow-code/DryCreek/DCEW-ICESat2-ATL06sr-atl08class-params';
ref_elevations_atl06_class = '/Users/karinazikan/Documents/ICESat2-snow-code/DryCreek/DCEW-ICESat2-ATL06sr-atl08class-ref-elevations';

%set colors
colors{1} = cmocean('-dense',5);
colors{2} = cmocean('-algae',4);
colors{3} = cmocean('ice',4);

%site abbreviation for file names
abbrev = 'DCEW';

%Set snowcover to 'snowon' or 'snowoff'
snowcover = 'snowoff';

%% Load data
%load the reference elevation data
% E_08 = readtable(ref_elevations_atl08);
E_06 = readtable(ref_elevations_atl06);
E_06_class = readtable(ref_elevations_atl06_class);

%load the ICESat-2 data
% I_08 = readtable(icesat2_atl08); %read in files
I_06 = readtable(icesat2_atl06);
I_06_class = readtable(icesat2_atl06_class);

footwidth = 11; % approx. width of icesat2 shot footprint in meters

% Plot Site + Tracks
DTM_name = '/Users/karinazikan/Documents/ICESat2-snow-code/DryCreek/DCEW-10mDEMcropped-tif/DCEW-DEMclip.tif';
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
x = x./(10^3); % to km
y = y./(10^3);
y = flip(y);% to km
DTM(DTM<0) = NaN;

fig1 = figure(1);
imagesc(x,y,DTM) 
daspect([1 1 1])
colormap(cmocean('grey'))
hold on
scatter([I_06.Easting./(10^3)],[I_06.Northing./(10^3)],[],'g','.')
xlabel('Easting (km)')
ylabel('Northing (km)')
set(gca,'fontsize',16);

% % Filter snow-on or snow-off for ATL08
% bright = I_08.Brightness_Flag;
% if snowcover == 'snowoff'
%     ib = find(bright == 0);
%     disp('Snow off')
% elseif snowcover == 'snowon'
%     ib = find(bright == 1);
%     disp('Snow on')
% else
%     error('snowcover must be set to snowon or snowoff')
% end
% I_08 = I_08(ib,:);
% E_08 = E_08(ib,:);

% Filter snow-on or snow-off for ATl06
dates_06 = datetime(I_06.time.Year,I_06.time.Month,I_06.time.Day);
if snowcover == 'snowoff'
    ib = find(I_06.time.Month >=6 & I_06.time.Month <= 9);
    disp('Snow off')
elseif snowcover == 'snowon'
    ib = find(I_06.time.Month <=2 & I_06.time.Month >=12);
    disp('Snow on')
else
    error('snowcover must be set to snowon or snowoff')
end
I_06 = I_06(ib,:);
E_06 = E_06(ib,:);

% Filter snow-on or snow-off for ATl06 atl08class
dates_06_class = datetime(I_06_class.time.Year,I_06_class.time.Month,I_06_class.time.Day);
if snowcover == 'snowoff'
    ib = find(I_06_class.time.Month >=6 & I_06_class.time.Month <= 9);
    disp('Snow off')
elseif snowcover == 'snowon'
    ib = find(I_06_class.time.Month <=2 & I_06_class.time.Month >=12);
    disp('Snow on')
else
    error('snowcover must be set to snowon or snowoff')
end
I_06_class = I_06_class(ib,:);
E_06_class = E_06_class(ib,:);


%% Loop through each product
N_Products = 3; %number of products to analyze
for i = 2:N_Products
    if i == 1
        T = I_08;
        E = E_08;
        acronym = 'ATL08';

        % ICESat-2 data
        zmod = T.Elevation_bestfit(:); % save the fitted 'model' elevations (icesat-2 elevations)
        %zstd = T.std; %save the standard deviation of the icesat-2 elevation estimates
        easts = T.Easting(:); % pull out the easting values
        norths = T.Northing(:); % pull out the northings
        slope = T.slope(:);
        aspect = T.aspect(:);
        %canopy = T.Canopy(:);
    elseif i == 2
         T = I_06;
        E = E_06;
        acronym = 'ATL06sr';    

        % ICESat-2 data
        zmod = T.h_mean(:); % save the median 'model' elevations (icesat-2 elevations
        zstd = T.h_sigma; %save the standard deviation of the icesat-2 elevation estimates
        easts = T.Easting(:); % pull out the easting values
        norths = T.Northing(:); % pull out the northings
        slope = T.slope(:);
        aspect = T.aspect(:);
%        canopy = T.Canopy(:);
    else
        T = I_06_class;
        E = E_06_class;
        acronym = 'ATL06sr atl08 classifications';
        % ICESat-2 data
        zmod = T.h_mean(:); % save the median 'model' elevations (icesat-2 elevations
        zstd = T.h_sigma; %save the standard deviation of the icesat-2 elevation estimates
        easts = T.Easting(:); % pull out the easting values
        norths = T.Northing(:); % pull out the northings
        slope = T.slope(:);
        aspect = T.aspect(:);
     %   canopy = T.Canopy(:);
    end
    % Ref elev data
    elevation_report(:,1) = E.elevation_report_nw_mean;  %non-weighted mean ref elevation
    elevation_report(:,2) = E.elevation_report_mean;    %weighted mean ref elevation
    elevation_report(:,3) = E.elevation_report_fitted;    %weighted & fitted ref elevation
    elevation_report_std = E.elevation_report_std;  %std of elevations within footprint

    %% Stats
    %calculate the elevation residuals
    for j = 1:3
        differences = zmod - elevation_report(:,j); %calculate the icesat2 elevations and the calculated reference elevations
        differences(differences > 80) = NaN; differences(differences < -80) = NaN; %remove extreme outliers
        Dmean(:,j) = nanmean(differences); % calculate mean of diferences
        Dstd(:,j) = std(differences,'omitnan'); % calculate std of diferences
        zrmse(:,j) = sqrt(nansum((differences).^2)./length(differences)); %calculate rmse of  differeces

        %         % Removing residuals below -13
        %         ix = find(differences < -13);
        %         zmod(ix) = NaN;
        %         differences(ix) = NaN;

        Residuals(:,j) = differences;
        %   Residuals(:,1) is comparison to non weighted mean elevations
        %   Residuals(:,2) is comparison to weighted mean elevations
        %   Residuals(:,3) is comparison to weighted & fitted elevations

    end
    Residuals(:,4) = i;
    ResidualsAll{i} = Residuals;
    %% Plots
    % Reference historgrams
    fig2 = figure(2);
    subplot(N_Products,1,i); set(gcf,'position',[50 50 800 500]); clear h;
    binwidth = 0.2;
    h(1) = histogram(Residuals(:,1),'Normalization','pdf'); h(1).BinWidth = binwidth; h(1).FaceAlpha = 1; h(1).FaceColor = colors{i}(1,:);  h(1).EdgeColor = 'k'; hold on;
    h(2) = histogram(Residuals(:,2),'Normalization','pdf'); h(2).BinWidth = binwidth; h(2).FaceAlpha = 0.75; h(2).FaceColor = colors{i}(2,:); h(2).EdgeColor = 'k';
    h(3) = histogram(Residuals(:,3),'Normalization','pdf'); h(3).BinWidth = binwidth; h(3).FaceAlpha = 0.75;  h(3).FaceColor = colors{i}(3,:); h(3).EdgeColor = 'w';
    plot([0,0],[0,.8], 'linewidth', 2, 'Color','k') % plot reference 0 line
    set(gca,'fontsize',16,'xlim',[-4 2]);
    legend(h,['DEM: Non-weighted mean, std = ' num2str(Dstd(:,1))],['DEM: Weighted mean, std = ' num2str(Dstd(:,2))],['DEM: Weighted and fitted, std = ' num2str(Dstd(:,3))],'Location','northwest');
    xlabel('Vertical offset (m)'); ylabel('Probability density'); title(acronym);

    % Reference pdfs
    fig3 = figure(3);
    subplot(N_Products,1,i); set(gcf,'position',[50 50 800 500]); hold on
    binwidth = 0.2;
    fplot(@(x) mynormpdf(x,Dmean(:,1), Dstd(:,1)),[-10 8], 'Linewidth', 2,'Color',colors{i}(1,:));
    fplot(@(x) mynormpdf(x,Dmean(:,2), Dstd(:,2)),[-10 8], 'Linewidth', 2,'Color',colors{i}(2,:));
    fplot(@(x) mynormpdf(x,Dmean(:,3), Dstd(:,3)),[-10 8], 'Linewidth', 2,'Color',colors{i}(3,:));
    plot([0,0],[0,.8], 'linewidth', 2, 'Color','k') % plot reference 0 line
    set(gca,'fontsize',16,'xlim',[-10 8]);
    legend(h,['DEM: Non-weighted mean, std = ' num2str(Dstd(:,1))],['DEM: Weighted mean, std = ' num2str(Dstd(:,2))],['DEM: Weighted and fitted, std = ' num2str(Dstd(:,3))],'Location','northwest');
    xlabel('Vertical offset (m)'); ylabel('Probability density'); title(acronym);

    % Reference pdfs for each tested product
    fig6 = figure(6); hold on
    if i == 1
        fplot(@(x) mynormpdf(x,Dmean(:,1), Dstd(:,1)),[-10 8], 'Linewidth', 2,'Color',colors{i}(3,:));
    else
        fplot(@(x) mynormpdf(x,Dmean(:,1), Dstd(:,1)),[-10 8], 'Linewidth', 2,'Color',colors{i}(3,:));
    end


    clear elevation_report Residuals
end
%% Plots outside loop
% Products historgrams
fig4 = figure(4); clf; hold on
for i = 2:N_Products
    if i == 1
        h(i) = histogram(ResidualsAll{i}(:,3),'Normalization','pdf');  h(i).FaceAlpha = 1; h(i).BinWidth = binwidth; h(i).FaceColor = colors{i}(3,:);  h(i).EdgeColor = 'k';
    else
        h(i) = histogram(ResidualsAll{i}(:,1),'Normalization','pdf');  h(i).FaceAlpha = .75; h(i).BinWidth = binwidth; h(i).FaceColor = colors{i}(3,:);  h(i).EdgeColor = 'k';
    end
end
plot([0,0],[0,.8], 'linewidth', 2, 'Color','k') % plot reference 0 line
set(gca,'fontsize',16,'xlim',[-4 2]);
set(gcf,'position',[50 50 800 400]);
legend('ATL08','ATL06sr','ATL06sr ATL08 clasifications','Location','northwest');
xlabel('Vertical offset (m)'); ylabel('Probability density');

% Finish reference pdfs for each tested product
fig6 = figure(6);
% plot([0,0],[0,.8], 'linewidth', 2, 'Color','k') % plot reference 0 line
set(gca,'fontsize',16);
set(gcf,'position',[50 50 800 400]);
legend('ATL08','ATL06sr','ATL06sr ATL08 clasifications','Location','northwest');
xlabel('Vertical offset (m)'); ylabel('Probability density');

%% create boxcharts for each terrain parameter
% ResidualsTable.residuals = [ResidualsAll{1}(:,3); ResidualsAll{2}(:,1); ResidualsAll{3}(:,1)];
% ResidualsTable.product = [ResidualsAll{1}(:,4); ResidualsAll{2}(:,4); ResidualsAll{3}(:,4)];
% ResidualsTable.elevations = [E_08.elevation_report_mean; E_06.elevation_report_mean; E_06_class.elevation_report_mean];
% ResidualsTable.aspect = [I_08.aspect; I_06.aspect; I_06_class.aspect];
% ResidualsTable.slope = [I_08.slope; I_06.slope; I_06_class.slope];

ResidualsTable.residuals = [ResidualsAll{2}(:,1); ResidualsAll{3}(:,1)];
ResidualsTable.product = [ResidualsAll{2}(:,4); ResidualsAll{3}(:,4)];
ResidualsTable.elevations = [E_06.elevation_report_mean; E_06_class.elevation_report_mean];
ResidualsTable.aspect = [I_06.aspect; I_06_class.aspect];
ResidualsTable.slope = [I_06.slope; I_06_class.slope];
figure(5);
%ELEVATION
h = histogram(E_08.elevation_report_mean(~isnan(ResidualsAll{1}(:,1))),5);
elev_binwidth = h.BinWidth; elev_binedges = h.BinEdges;
clear h;
% close(gcf);
%SLOPE
h = histogram(I_08.slope(~isnan(ResidualsAll{1}(:,1))),5);
slope_binwidth = h.BinWidth; slope_binedges = h.BinEdges;
clear h;
% close(gcf);
%ASPECT
h = histogram(I_08.aspect(~isnan(ResidualsAll{1}(:,1))),5);
aspect_binwidth = h.BinWidth; aspect_binedges = h.BinEdges;
clear h;
% close(gcf);

whiskerline = '-'; outliermarker = 'o';

% boxplot figure
fig5 = figure(5); clf
%ELEVATION
subplot(3,1,1);
hold on
bins = {num2str(elev_binedges(2)) num2str(elev_binedges(3)) num2str(elev_binedges(4)) num2str(elev_binedges(5)) num2str(elev_binedges(6))};
groupElev = discretize(ResidualsTable.elevations,elev_binedges,'categorical',bins);
%meanWeight = groupsummary(ResidualsTable.residuals,groupElev,'mean');
boxchart(groupElev,ResidualsTable.residuals,'GroupByColor',ResidualsTable.product,'MarkerStyle','none')
colororder([colors{1}(3,:); colors{2}(3,:); colors{3}(3,:)]); 
ylim([-5,3])
% plot(meanWeight,'-o','Color','k')
set(gca,'fontsize',16,'box','on'); drawnow;
legend('ATL08','ATL06sr','ATL06sr ATL08 clasifications','Location','northwest');
xlabel('Elevation (m a.s.l.)','fontsize',16); %ylabel('Elevation residuals (m)','fontsize',16);
%text(1,max(ylims)-0.05*range(ylims),'a)','fontsize',16);
%ASPECT
subplot(3,1,2);
hold on
bins = {num2str(aspect_binedges(2)) num2str(aspect_binedges(3)) num2str(aspect_binedges(4)) num2str(aspect_binedges(5)) num2str(aspect_binedges(6))};
groupAspect = discretize(ResidualsTable.aspect,aspect_binedges,'categorical',bins);
%meanWeight = groupsummary(ResidualsTable.residuals,groupAspect,'mean');
boxchart(groupAspect,ResidualsTable.residuals,'GroupByColor',ResidualsTable.product,'MarkerStyle','none')
ylim([-5,3])
% plot(meanWeight,'-o','Color','k')
set(gca,'fontsize',16,'box','on'); drawnow;
xlabel('Aspect (degrees)','fontsize',16); ylabel('Elevation residuals (m)','fontsize',16);
%text(1,max(ylims)-0.05*range(ylims),'a)','fontsize',16);
%SLOPE
subplot(3,1,3);
hold on
bins = {num2str(slope_binedges(2)) num2str(slope_binedges(3)) num2str(slope_binedges(4)) num2str(slope_binedges(5)) num2str(slope_binedges(6))};
groupSlope = discretize(ResidualsTable.aspect,slope_binedges,'categorical',bins);
%meanWeight = groupsummary(ResidualsTable.residuals,groupAspect,'mean');
boxchart(groupSlope,ResidualsTable.residuals,'GroupByColor',ResidualsTable.product,'MarkerStyle','none')
ylim([-5,3])
% plot(meanWeight,'-o','Color','k')
set(gca,'fontsize',16,'box','on'); drawnow;
xlabel('Slope (degrees)','fontsize',16); %ylabel('Elevation residuals (m)','fontsize',16);
%text(1,max(ylims)-0.05*range(ylims),'a)','fontsize',16);


