%%% This code was writen by Karina for the GEOPH 522 final project and is built off
%%% pervious codes writen by Colten and Ellyn, batch_icesat2_coregistration.m and coregister_icesat2
%%% The inputs and initial window calcultions were coded by Colten, the
%%% start of Karina's code is marked by the line "Karina's Code Starts HERE"


%% Inputs
clearvars; close all;
addpath(['./functions']) 

%DTM (be sure the path ends in a /)
DTM_path = 'RCEW_DEM/';
DTM_name = 'RCEW_1m_WGS84UTM11_WGS84.tif';
if contains(DTM_name,'.tif')
    DTM_date = '20120826'; %only need to change this if the DTM is a geotiff
end

%csv (be sure the path ends in a /)
csv_path = '/Users/karinazikan/Documents/ICESat2-snow-code/';

%site abbreviation for file names
abbrev = 'RCEW';

%ICESat-2 product acronym
acronym = 'ATL08';
ATL0X = 8;

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
% for i = 1:length(csvs)
%     icesat2 = [csv_path,'RCEW-ICESat2-ATL08-params']; %compile the file name
%     file = readtable(icesat2); %read in files
%     T = [T; file]; %combine tables
% end
icesat2 = [csv_path,'RCEW-ICESat2-ATL08-params']; %compile the file name
file = readtable(icesat2); %read in files
T = [T; file];

T = T([1:250],:);
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

% for ATL08 files only use snow-free data (brightness flag == 0)
if contains(icesat2, 'ATL08') % ATL08 commands
    bright = T.Brightness_Flag;
    ib = find(bright == 0);
    easts = easts(ib);
    norths = norths(ib);
    zmod = zmod(ib);
    zmodfit = zmodfit(ib);
end

cd '/Users/karinazikan/Documents/ICESat2-snow-code'

%% Calculating footprints and reference elevation for each data point
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
%%
for r=1:length(zmod)
    
    %identify the R2erence elevation points in each ICESat2 footprint
    xv = xc(r,[3:6 3]); % bounding box x vector
    yv = yc(r,[3:6 3]); % bounding box y vector

    %data in the footprint
    in = inpolygon(xgrid, ygrid, xv, yv); % get logical array of in values
    pointsinx = xgrid(in); % save x locations
    pointsiny = ygrid(in); % save y locations
    elevationsin = elevations(in); % save elevations
    %elevation_report(t) = nanmean(elevationsin);
    

    %wieghted average
    dist = nan([1,length(pointsinx)])'; %initialize dist
    for a = 1:length(pointsinx)
        phi = atan2d((pointsiny(a)-norths(r)),(pointsinx(a)-easts(r)));
        if phi > theta
            dist(a)=abs(sqrt((pointsiny(a)-norths(r))^2+(pointsinx(a)-easts(r))^2)*sind(phi-theta(r))); %distance from the line in the center of the window
        else
            dist(a)=abs(sqrt((pointsiny(a)-norths(r))^2+(pointsinx(a)-easts(r))^2)*sind(phi+theta(r)));
        end
    end

    maxdist = footwidth/2; % defining the maximum distance a point can be from the center icesat2 point
    w = 15/16*(1-(dist/maxdist).^2).^2; %bisqared kernel
    elevation_report(r,:) = sum(w.*elevationsin)./sum(w); %unbiased estimate
    %     elevation_report(r,:) = mean(elevationsin);
    elevation_report_std(r,:) = std(elevationsin); %std of the elevations within the footprint

    %non wieghted average
    elevation_report_nw(r,:) = nanmean(elevationsin); % non-wieghted elevations

    %weighted fit
    warning('off')
    p{1} = fit([pointsinx, pointsiny],elevationsin,'poly11','Weights',w); %fit linear polynomial
    p{2} = fit([pointsinx, pointsiny],elevationsin,'poly33','Weights',w); %fit cubic polynomial
    p{3} = fit([pointsinx, pointsiny],elevationsin,'poly44','Weights',w); %fit quadratic polynomial
    %     p{1} = fit([pointsinx, pointsiny],elevationsin,'poly11'); %fit linear polynomial
    %     p{2} = fit([pointsinx, pointsiny],elevationsin,'poly33'); %fit cubic polynomial
    %     p{3} = fit([pointsinx, pointsiny],elevationsin,'poly44'); %fit quadratic polynomial
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
    elevation_report_fitted(r,:) = Em(a);
end
%% Stats

%calculate the elevation residuals
differences = zmod - elevation_report; %calculate the difference between the mean icesat2 elevations and the calculated reference elevations
differences(differences > 80) = nan; differences(differences < -80) = nan; %remove extreme outliers
Dmean = nanmean(differences); % calculate mean of diferences
Dstd = std(differences,'omitnan'); % calculate std of diferences
zrmse = sqrt(nansum((differences).^2)./length(differences)); %calculate rmse of  differeces

differencesfit = zmodfit - elevation_report_fitted; %calculate the difference between the fitted icesat2 elevations and the calculated reference elevations
differencesfit(differencesfit > 80) = nan; differencesfit(differencesfit < -80) = nan; %remove extreme outliers
Dfitmean = nanmean(differencesfit); % calculate mean of fit diferences
Dfitstd = std(differencesfit,'omitnan'); % calculate std of diferences
zfitrmse = sqrt(nansum((differencesfit).^2)./length(differencesfit)); %calculate rmse of fit differeces

% Removing residuals below -13
ix = find(differences < -13);
zmod(ix) = NaN;
differences(ix) = NaN;
ix = find(differencesfit < -13);
zmodfit(ix) = NaN;
differencesfit(ix) = NaN;
%% Plots
% 1-1 plots
fig1 = figure(1); clf %create figure 1
hold on
errorbar(elevation_report,zmod,zstd,zstd,elevation_report_std,elevation_report_std,'.','MarkerSize',10); % 1-1 plot of median icesat2 elev vs reference elev
plot([1000 2100], [1000 2100], '-r') % plot reference 1-1 line
xlabel('reference elevation'); %label x axis
ylabel('ICESat-2 mean elevation') %label y axis
legend('Data', 'Reference Line') %legendhold off

text(1100,2100,['RMSE = ' num2str(zrmse)]) %print rmse on plot

fig2 = figure(2); clf %create figure 2
hold on
errorbar(elevation_report_fitted,zmodfit,zstd,zstd,elevation_report_std,elevation_report_std,'.','MarkerSize',10); % 1-1 plot of fitted icesat2 elev vs reference elev
xlabel('reference elevation'); %label x axis
ylabel('ICESat-2 fitted elevation') %label y axis
plot([1000 2100], [1000 2100], '-r') % plot reference 1-1 line
hold off
legend('Data', 'Reference Line') %legend
text(1100,2100,['RMSE = ' num2str(zfitrmse)]) %print rmse on plot

% Histogram of diferences
nbins = 30; %sets the number of bins for the reletive dencity histograms (This value will also be used for the plots in fig 2)

% plot mean diferences
fig3 = figure(3); clf % open figure 3
[N,binx] = myRelDencHist(differences,nbins); %Calculates a relitive dencity histogram of the differences between the mean elevation and the DEM
relhist = bar(binx,N,1); %plots the relative dencity histogram for diferences
xlabel('Difference between the mean ICESat_2 elevation and the reference DEM'); % labeling the x axis
ylabel('Normaized number of observations'); % labeling the y axis
hold on %puts hold on so the histogram is not over written by the pdf plot
fplot(@(x) mynormpdf(x,Dmean, Dstd),[binx(1,1) binx(1,30)], 'Linewidth', 2); %uses the "true" mean and std for elevations to plot the normal pdf for elevations
plot([0,0],[0,max(N)+.05], 'linewidth', 2) % plot reference 0 line
hold off %turns off hold

% plot fitted diferences
fig4 = figure(4); clf % open figure 4
[N,binx] = myRelDencHist(differencesfit,nbins); %Calculates a relitive dencity histogram of the differences between the mean elevation and the DEM
relhist = bar(binx,N,1); %plots the relative dencity histogram for diferences
xlabel('Difference between the fitted ICESat_2 elevation and the reference DEM'); % labeling the x axis
ylabel('Normaized number of observations'); % labeling the y axis
hold on %puts hold on so the histogram is not over written by the pdf plot
fplot(@(x) mynormpdf(x,Dfitmean, Dfitstd),[binx(1,1) binx(1,30)], 'Linewidth', 2); %uses the "true" mean and std for elevations to plot the normal pdf for elevations
plot([0,0],[0,max(N)+.05], 'linewidth', 2) % plot reference 0 line
hold off %turns off hold

fig10 = figure(10); clf % open figure 3
[N,binx] = myRelDencHist(differences,nbins); %Calculates a relitive dencity histogram of the differences between the mean elevation and the DEM
relhist = bar(binx,N,1,'FaceAlpha',.75); %plots the relative dencity histogram for diferences
xlabel('Elevation Residuals'); % labeling the x axis
ylabel('Normaized number of observations'); % labeling the y axis
hold on %puts hold on so the histogram is not over written by the pdf plot
[N,binx] = myRelDencHist(differencesfit,nbins); %Calculates a relitive dencity histogram of the differences between the mean elevation and the DEM
relhist = bar(binx,N,1,'FaceAlpha',.75); %plots the relative dencity histogram for diferences
legend('Mean Elevation', 'Fitted Elevation')
set(gca,'FontSize',16)

% boxplot of differences
fig5 = figure(5); clf % open and clear fig 5
boxplot([differences,differencesfit],'Notch','on','Labels',{'mean','fitted'}) % boxplot
set(gca,'FontSize',16)
% 
% % plot of ICESat2 elevations and DEM vs norths
% fig6 = figure(6); clf %open and clear fig 6
% plot(norths([1:251],:),elevation_report([1:251],:),norths([1:251],:),zmod([1:251],:),'LineWidth',2); %plot elevation track of reference elevations and icesat2
% xlabel('Northing'); % labeling the x axis
% ylabel('Elevation'); % labeling the y axis

% DEM fit Examples
figure(7);
scatter3(pointsinx,pointsiny,elevationsin,'.','m')
hold on
plot(p{1})
xlabel('Easting (km)')
ylabel('Northing (km)')
zlabel('Elevation (m)')
colormap gray
set(gca,'FontSize',16)

figure(8);
scatter3(pointsinx,pointsiny,elevationsin,'.','m')
hold on
plot(p{2})
xlabel('Easting (km)')
ylabel('Northing (km)')
zlabel('Elevation (m)')
colormap gray
set(gca,'FontSize',16)

figure(9);
scatter3(pointsinx,pointsiny,elevationsin,'.','m')
hold on
plot(p{3})
xlabel('Easting (km)')
ylabel('Northing (km)')
zlabel('Elevation (m)')
colormap gray
set(gca,'FontSize',16)

% %% Uncertanty in mean error
% nMC = 1000; %defines the number of times to do the Monte Carlo simulation
% nsamp = 5:1:70; %makes a list of sample sizes to do the Monte Carlo for
% 
% for n = 1:length(nsamp) %creates for loop to varry sample size
%     for n2 = 1:nMC %creates a for loop to do monte carlo simulations
%         D2 = randsample(differences, nsamp(n), true); %Samples 10 random datapoints from the differences
%         D_std(n2) = std(D2); %Stores standard deviation of D2
%         Dmu(n2) = mean(D2); %Stores mean of D2
%         D3 = randsample(differencesfit, nsamp(n), true); %Samples 10 random datapoints from the differences
%         Dfit_std(n2) = std(D3); %Store standard deviation of D3
%         Dfitmu(n2) = mean(D3); %Stores mean of D3
%     end %ends loop
%     Dstd_mu(n) = mean(Dstd); %stores the mean sample stds for a sample of size nsamp(n)
%     Dstd_std(n) = std(Dstd); %stores the std sample stds for a sample of size nsamp(n)
%     Dmu_mu(n) = mean(Dmu); %stores the mean sample means for a sample of size nsamp(n)
%     Dmu_std(n) = std(Dmu); %stores the std sample means for a sample of size nsamp(n)
%     Dfitstd_mu(n) = mean(Dfitstd); %stores the mean sample stds for a sample of size nsamp(n)
%     Dfitstd_std(n) = std(Dfitstd); %stores the std sample stds for a sample of size nsamp(n)
%     Dfitmu_mu(n) = mean(Dfitmu); %stores the mean sample means for a sample of size nsamp(n)
%     Dfitmu_std(n) = std(Dfitmu); %stores the std sample means for a sample of size nsamp(n)
% end %ends loop
% 
% %% Varying sample size plot
% % mean icesat2 elevations
% fig7 = figure(7); clf %creates and clears figure 7
% 
% subplot(1,2,1); %subplot 1
% plot(nsamp, Dstd_mu,nsamp,Dstd_mu + Dstd_std,nsamp,Dstd_mu - Dstd_std); %plots mean of sample std for different sample sizes and mean min plus/minus std
% hold on;  %puts hold on
% legend('mean difference std','mean difference std + std', 'mean difference std - std') %creates legend
% xlabel('Sample size'); % labeling the x axis
% ylabel('Mean std of differeces between icesat2 mean elevations and the reference DEM'); % labeling the y axis
% 
% subplot(1,2,2); %subplot 2
% plot(nsamp, Dmu_mu,nsamp,Dmu_mu + Dmu_std,nsamp,Dmu_mu - Dmu_std); %plots mean of sample mean for different sample sizes and mean min plus/minus std
% hold on;  %puts hold on
% legend('mean difference mean','mean difference mean + std', 'mean difference mean - std') %creates legend
% xlabel('Sample size'); % labeling the x axis
% ylabel('Mean mean of differeces between icesat2 mean elevations and the reference DEM'); % labeling the y axis
% 
% % fitted icesat2 elevations
% fig8 = figure(8); clf %creates and clears figure 8
% 
% subplot(1,2,1); %subplot 1
% plot(nsamp, Dfitstd_mu,nsamp,Dfitstd_mu + Dfitstd_std,nsamp,Dfitstd_mu - Dfitstd_std); %plots mean of sample std for different sample sizes and mean min plus/minus std
% hold on;  %puts hold on
% legend('mean difference std','mean difference std + std', 'mean difference std - std') %creates legend
% xlabel('Sample size'); % labeling the x axis
% ylabel('Mean std of differeces between icesat2 fitted elevations and the reference DEM'); % labeling the y axis
% 
% subplot(1,2,2); %subplot 2
% plot(nsamp, Dfitmu_mu,nsamp,Dfitmu_mu + Dfitmu_std,nsamp,Dfitmu_mu - Dfitmu_std); %plots mean of sample mean for different sample sizes and mean min plus/minus std
% hold on;  %puts hold on
% legend('mean difference mean','mean difference mean + std', 'mean difference mean - std') %creates legend
% xlabel('Sample size'); % labeling the x axis
% ylabel('Mean mean of differeces between icesat2 fitted elevations and the reference DEM'); % labeling the y axis

%% Saving Figures
% print(fig1,'-djpeg','Mean_1to1.jpeg'); %saves figure 1
% print(fig2,'-djpeg','Fitter_1to1.jpeg'); %saves figure 2
% print(fig3,'-djpeg','Mean_DifHist.jpeg'); %saves figure 3
% print(fig4,'-djpeg','Fitted_DifHist.jpeg'); %saves figure 4
% print(fig5,'-djpeg','boxplots.jpeg'); %saves figure 5
% print(fig6,'-djpeg','Transect.jpeg'); %saves figure 6
% print(fig7,'-djpeg','Mean_DifUncercertanty.jpeg'); %saves figure 7
% print(fig8,'-djpeg','Fitted_DifUncercertanty.jpeg'); %saves figure 8
