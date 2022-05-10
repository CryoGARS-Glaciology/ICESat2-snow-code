%%% This code was writen by Karina Zikan to statistically compare ICESat-2
%%% elevations compiled by batch_icesat2_terrain_parameter_comparison.m to 
%%% the reference elevations calculated by DTM_reference_elevations_calculation.m
%%%
%%% SPECIFIED INPUTS:
%%%     icesat2 = path to ICESat-2 file (output from batch_icesat2_terrain_parameter_comparison.m) 
%%%     ref_elevations = path to reference elevations file (output from DTM_reference_elevations_calculation.m)
%%%     abbrev = site abriviation for file name
%%%     acronym = ICESat-2 product acronym
%%%     snowcover = set to snow-on ('snowon') or snow-off ('snowoff') conditions
%%% OUTPUTS:
%%%     
%%%
%%% Last updated: May 10 2022 by Karina Zikan

%% Inputs
clearvars; close all;
addpath(['./functions']) 

%File paths
icesat2 = '/Users/karinazikan/Documents/ICESat2-snow-code/RCEW-ICESat2-ATL08-params';
ref_elevations = '/Users/karinazikan/Documents/ICESat2-snow-code/RCEW-ICESat2-ATL08-ref-elevations';

%site abbreviation for file names
abbrev = 'RCEW';

%ICESat-2 product acronym
acronym = 'ATL08';

%Set snowcover to 'snowon' or 'snowoff'
snowcover = 'snowoff';

%% Load data
%load the reference elevation data
E = readtable(ref_elevations);

%load the ICESat-2 data
T  = readtable(icesat2); %read in files
footwidth = 11; % approx. width of icesat2 shot footprint in meters

% Filter snow-on or snow-off for ATL08
if acronym == 'ATL08' % ATL08 commands
    bright = T.Brightness_Flag;
    if snowcover == 'snowoff'; 
        ib = find(bright == 0);
        disp('Snow off')
    elseif snowcover = ='snowon';
        ib = find(bright == 1);
    else
        error('snowcover must be set to snowon or snowoff')
        disp('Snow on')
    end
    T = T(ib);
    E = E(ib);
end

% ICESat-2 data  
zmod = T.Elevation(:); % save the median 'model' elevations (icesat-2 elevations)
zmodfit = T.Elevation_bestfit(:); % save the fitted 'model' elevations (icesat-2 elevations_bestfit)
zmodfit(isnan(zmod)) = NaN;
zstd = T.std; %save the standard deviation of the icesat-2 elevation estimates
easts = T.Easting(:); % pull out the easting values
norths = T.Northing(:); % pull out the northings
% Ref elev data
elevation_report_nw_mean = E.elevation_report_nw_mean;
elevation_report_mean = E.elevation_report_nw_mean;
elevation_report_fitted = E.elevation_report_fitted;
elevation_report_std = E.elevation_report_std;

%% Stats
%calculate the elevation residuals
differences = zmod - elevation_report_mean; %calculate the difference between the mean icesat2 elevations and the calculated reference elevations
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
errorbar(elevation_report_mean,zmod,zstd,zstd,elevation_report_std,elevation_report_std,'.','MarkerSize',10); % 1-1 plot of median icesat2 elev vs reference elev
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
% plot(norths([1:251],:),elevation_report_mean([1:251],:),norths([1:251],:),zmod([1:251],:),'LineWidth',2); %plot elevation track of reference elevations and icesat2
% xlabel('Northing'); % labeling the x axis
% ylabel('Elevation'); % labeling the y axis