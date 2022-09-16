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

%File paths
icesat2_atl08 = '/Users/karinazikan/Documents/ICESat2-snow-code/RCEW/RCEW-ICESat2-ATL08-params';
ref_elevations_atl08 = '/Users/karinazikan/Documents/ICESat2-snow-code/RCEW/RCEW-ICESat2-ATL08-ref-elevations';

icesat2_atl06_class = '/Users/karinazikan/Documents/ICESat2-snow-code/RCEW/RCEW-ICESat2-ATL06sr-atl08class';
ref_elevations_atl06_class = '/Users/karinazikan/Documents/ICESat2-snow-code/RCEW/RCEW-ICESat2-ATL06sr-atl08class-ref-elevations';

%site abbreviation for file names
abbrev = 'RCEW';

%Set snowcover to 'snowon' or 'snowoff'
snowcover = 'snowoff';

%% Load data
%load the reference elevation data
E_08 = readtable(ref_elevations_atl08);
E_06_class = readtable(ref_elevations_atl06_class);

%load the ICESat-2 data
I_08  = readtable(icesat2_atl08); %read in files
I_06_class = readtable(icesat2_atl06_class);

footwidth = 11; % approx. width of icesat2 shot footprint in meters

% Filter snow-on or snow-off for ATL08
bright = I_08.Brightness_Flag;
if snowcover == 'snowoff';
    ib = find(bright == 0);
    disp('Snow off')
elseif snowcover == 'snowon';
    ib = find(bright == 1);
    disp('Snow on')
else
    error('snowcover must be set to snowon or snowoff')
end
I_08 = I_08(ib,:);
E_08 = E_08(ib,:);

% Filter snow-on or snow-off for ATl06
dates_06_class = datetime(I_06_class.time.Year,I_06_class.time.Month,I_06_class.time.Day);
if snowcover == 'snowoff';
    ib = find(I_06_class.time.Month >=6 & I_06_class.time.Month <= 9);
    disp('Snow off')
elseif snowcover == 'snowon';
    ib = find(I_06_class.time.Month <=2 & I_06_class.time.Month >=12);
    disp('Snow on')
else
    error('snowcover must be set to snowon or snowoff')
end
I_06_class = I_06_class(ib,:);
E_06_class = E_06_class(ib,:);


%% Loop through each product
N_Products = 2; %number of products to analyze 
for i = 1:N_Products
    if i == 1
        T = I_08;
        E = E_08;
        acronym = 'ATL08';

        % ICESat-2 data
        zmod = T.Elevation_bestfit(:); % save the fitted 'model' elevations (icesat-2 elevations)
        zstd = T.std; %save the standard deviation of the icesat-2 elevation estimates
        easts = T.Easting(:); % pull out the easting values
        norths = T.Northing(:); % pull out the northings
        slope = T.slope(:);
        aspect = T.aspect(:);
        canopy = T.Canopy(:);
    else %i == 2
        T = I_06_class;
        E = E_06_class;
        acronym = 'ATL06sr_class';

        % ICESat-2 data
        zmod = T.h_mean(:); % save the median 'model' elevations (icesat-2 elevations
        zstd = T.h_sigma; %save the standard deviation of the icesat-2 elevation estimates
        easts = T.Easting(:); % pull out the easting values
        norths = T.Northing(:); % pull out the northings
        %     slope = T.slope(:);
        %     aspect = T.aspect(:);
        %     canopy = T.Canopy(:);

        %     else
        %         T = I_06;
        %         E = E_06;
        %         acronym = 'ATL06sr';
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

    %% Plots
    % Reference historgrams
    figure(1); 
    subplot(N_Products,1,i); clf; set(gcf,'position',[50 50 800 400]); clear h;
    binwidth = 0.2;
    h(1) = histogram(Residuals(:,1),'Normalization','pdf'); h(1).BinWidth = binwidth; h(1).FaceAlpha = 1;  h(1).EdgeColor = 'k'; hold on;
    h(2) = histogram(Residuals(:,2),'Normalization','pdf'); h(2).BinWidth = binwidth; h(2).FaceAlpha = 0.75;  h(2).EdgeColor = 'k'; 
    h(3) = histogram(Residuals(:,3),'Normalization','pdf'); h(3).BinWidth = binwidth; h(3).FaceAlpha = 0.75;  h(3).EdgeColor = 'w'; 
    plot([0,0],[0,.7], 'linewidth', 2, 'Color','k') % plot reference 0 line
    set(gca,'fontsize',16,'xlim',[nanmedian(T.VerticalErrors)-3*1.4826*mad(T.VerticalErrors,1) nanmedian(T.VerticalErrors)+3*1.4826*mad(T.VerticalErrors,1)]);
    leg = legend(h,'non-weighted mean','weighted mean','weighted & fitted');
    xlabel('Vertical offset (m)'); ylabel('Probability density');

end