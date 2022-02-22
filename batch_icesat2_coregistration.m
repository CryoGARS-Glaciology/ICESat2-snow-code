%%%
%Loop through all the csv files containing ICESat-2 ATL08 or ATL06 data and
%coregister to a reference snow-off digital terrain model (DTM). You must
%specify the following in the section below:
%1) Path where you keep your codes 
%2) Path for the DTM
%3) Name of the DTM (if analyzing elevations for a glacier, create a
%timeseries of glacier elevations using concatenate_glacier_DEM_timeseries.m)
%4) Path for the csv files
%5) Site abbreviation used for compiled outputs
%6) ICESat-2 data acronym (e.g., ATL06, ATL08)
%7) Path with filename for glacier outline (if working with ATL06 data)

%Outputs: 
%1) A matrix, Abest, containing the easting and northing offsets that yield
%the smallest root mean square difference in elevation.
%2) A matrix, RMSDbest, containing the minimum root mean squared difference in
%elevations between the DTM and ICESat-2 data yielded by the horizontal
%coregistration process.
%3) Revised csv files with '-edited' appended to the raw csv filenames that
%contain the horizontally & vertically coregistered data.
%4) A csv table, [abbrev,'-ICESat2-',acronym,'-params.csv'], with all the
%coregistered transects and season flags, coregistered ICESat-2 elevations,
%and differences between the coregistered ICESat-2 and reference ground elevations.

%%%
%% Inputs
clearvars; close all; 

%path for the local version of the Github directory
addpath('/Users/ellynenderlin/Research/NASA_CryoIdaho/ICESat2-snow-code/')

%DTM (be sure the path ends in a /)
DTM_path = '/Users/ellynenderlin/Research/NASA_CryoIdaho/mountains/RCEW/DEMs/';
DTM_name = 'RCEW_1m_WGS84UTM11_WGS84.tif';
if contains(DTM_name,'.tif')
    DTM_date = '20071114'; %only need to enter datestring in file name if the reference elevation map is a geotiff
end

%ROI polygon in UTM coordinates (not necessary for ATL08 data)
S = shaperead('/Users/ellynenderlin/Research/NASA_CryoIdaho/mountains/RCEW/ROIs/RCEW-outline_WGS84-UTM11N.shp'); %glacier outline if using ATL06

%csv (be sure the path ends in a /)
csv_path = '/Users/ellynenderlin/Research/NASA_CryoIdaho/mountains/RCEW/csvs/';

%site abbreviation for file names
abbrev = 'RCEW'; 

%ICESat-2 product acronym
acronym = 'ATL08';

%if ATL06 equilbrium line altitude or typical late summer snowline
%if ATL08 typical rain-snow transition line in catchment (need for coregistration if not enough summer data)
if contains(acronym,'ATL06') %SITE SPECIFIC!!!
   snowline = 1235; %Wolverine Glacier = 1235; 
else %typical rain-snow transition elevation (~1500 m for contiguous US Rockies)
   snowline = 1400; %RCEW = 1400;
end

%days of year
modays_norm = [31 28 31 30 31 30 31 31 30 31 30 31];
cumdays_norm = cumsum(modays_norm); cumdays_norm = [0 cumdays_norm(1:11)];
modays_leap = [31 29 31 30 31 30 31 31 30 31 30 31];
cumdays_leap = cumsum(modays_leap); cumdays_leap = [0 cumdays_leap(1:11)];

%% Coregistration (ONLY RUN ONCE OR YOU NEED TO MODIFY THE CSV NAME SEARCH ON LINE 69)
%read in the snow-off reference elevation map
cd_to_DTM = ['cd ',DTM_path]; eval(cd_to_DTM);
disp('Loading DTM(s)');
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

%loop through the csvs & determine the horizontal offsets
tic;
for i = 1:length(csvs)
    disp(['Coregistering csv #',num2str(i),' (',csvs(i).name,')']);
    ICESat2_file = [csv_path,csvs(i).name]; %compile the file name
    
    %if analyzing ATL06 data for a glacier, select the most recent DEM
    if contains(acronym,'ATL06')
%         clear DTM Ref;

        %convert datestamp to decimal dates
        yr = str2double(csvs(i).name(7:10));
        if mod(yr,4) == 0
            decidate = (cumdays_leap(str2double(csvs(i).name(11:12)))+str2double(csvs(i).name(13:14)))/sum(modays_leap);
        else
            decidate = (cumdays_norm(str2double(csvs(i).name(11:12)))+str2double(csvs(i).name(13:14)))/sum(modays_norm);
        end
        deciyear = yr+decidate;
        
        %identify the most recent DEM
        if contains(DTM_name,'.mat')
            yr_ref = find(DEMdate<=deciyear,1,'last');
            DTM = Z(yr_ref).z; Ref = Z(yr_ref).R;
            disp(['Using DEM from ',num2str(Z(yr_ref).deciyear)]);
            clear yr_ref;
        end
        clear yr decidate deciyear;
    end
    
    %determine the map-view offset
    myFuncHandle = @(A)coregister_icesat2(ICESat2_file,DTM,Ref,A); %create the handle to call the coregistration function
    vals = readmatrix(ICESat2_file); %read in the csv to check its size
    if size(vals,1) > 1
        [Abest(i,:),RMSDbest(i)] = fminsearch(myFuncHandle,[0,0]); %initial horizontal offset estimate = [0,0] = [0 m East, 0 m North]
    else
        Abest(i,:) = [NaN NaN]; RMSDbest(i) = NaN;
    end
    save('Abest.mat','Abest','-v7.3'); save('RMSDbest.mat','RMSDbest','-v7.3');
    disp(['x-offset = ',num2str(Abest(i,1)),' m & y-offset = ',num2str(Abest(i,2)),' m w/ RMSD = ',num2str(RMSDbest(i)),' m']);
    clear vals;
    toc; tic; %display processing time
    
    %horizontally coregister & write the reference DTM and vertical offset
    %values to the csv file
    T = readtable(csvs(i).name);
    %specify variables with column indices in square brackets if you want to cull data
    if contains(csvs(i).name, 'ATL08')
        t = T(:,:); 
    else
        t = T(:,:);
    end
    %apply horizontal coregistration shift that minimizes the rmsd in elevations determined using fminsearch above
    if ~isnan(Abest(i,1))
        t.Northing = t.Northing+Abest(i,2); t.Easting = t.Easting+Abest(i,1); 
    end
    %extract reference elevations and vertical errors
    %NOTE: vertical errors = ICESat2-reference so errors > 0 mean ICESat2
    %elevations are greater than the reference surface (possible snow signal)
    if length(t.Elevation) > 1
        [t.ReferenceElevation,t.VerticalErrors,~] = extract_icesat2_vertical_errors(ICESat2_file,DTM,Ref);
        disp(['Median vertical bias after horizontal coregistration = ',num2str(nanmedian(t.VerticalErrors)),' m']);
    else
        t.ReferenceElevation = NaN; t.VerticalErrors = NaN;
        disp('Too few data to coregister!');
    end
    writetable(t,[csvs(i).name(1:end-4),'-edited.csv']);
    disp(['Adjusted coordinates and reference elevations saved to ',csvs(i).name(1:end-4),'-edited.csv']);
    clear T t;
    disp(' '); %leave a space in the command line after the outputs
end

%% apply median snow-free offset to vertically coregister 

%identify the ICESat-2 data csv files
cd_to_csv = ['cd ',csv_path]; eval(cd_to_csv);
csvs = dir('*-edited.csv');

%days of year for calculation of decimal dates
modays_norm = [31 28 31 30 31 30 31 31 30 31 30 31];
cumdays_norm = cumsum(modays_norm); cumdays_norm = [0 cumdays_norm(1:11)];
modays_leap = [31 29 31 30 31 30 31 31 30 31 30 31];
cumdays_leap = cumsum(modays_leap); cumdays_leap = [0 cumdays_leap(1:11)];

%loop through the csvs, assign seasons, & compile into one table
seasonal_bias = []; seasonal_id = []; seasonal_refz = [];
disp('Extracting vertical offsets from snow-off data...');
for i = 1:length(csvs)
    ICESat2 = [csv_path,csvs(i).name]; %compile the file name
    t = readtable(csvs(i).name);

    %identify season (1=winter(Jan,Feb,Mar),2=spring(Apr,May,Jun),3=summer(Jul,Aug,Sep),4=fall(Oct,Nov,Dec))
    seas = csvs(i).name(11:12); 
    if str2double(seas) >= 1 && str2double(seas) <= 3
        seas_flag = 1;
    elseif str2double(seas) >= 4 && str2double(seas) <= 6
        seas_flag = 2;
    elseif str2double(seas) >= 7 && str2double(seas) <= 9
        seas_flag = 3;
    else
        seas_flag = 4;
    end
    t.season = repmat(seas_flag,size(t.Elevation));
    
    %compile seasonal elevation biases
    if contains(csvs(i).name, 'ATL08')
%         if str2double(seas) >= 6 && str2double(seas) <= 11 %based on SNOTEL snow-free months
            seasonal_bias = [seasonal_bias; t.VerticalErrors(t.Brightness_Flag==0)]; %make sure only snow-free elevation differences are used
            seasonal_id = [seasonal_id; t.season(t.Brightness_Flag==0)];
            seasonal_refz = [seasonal_refz; t.ReferenceElevation(t.Brightness_Flag==0)];
%         end
    else
%         if str2double(seas) >= 7 && str2double(seas) <= 9 %based on typical Alaska seasonality
            vert_bias = t.VerticalErrors(t.ReferenceElevation<=snowline); %must exclude data above snowline year-round for the glacier since the snowline marks the elevation of year-round snow
            vert_bias_season = t.season(t.ReferenceElevation<=snowline); vert_bias_refz = t.ReferenceElevation(t.ReferenceElevation<=snowline);
            in = inpolygon(t.Easting, t.Northing, S.X, S.Y); %identify footprints in the glacier outline
            vert_bias(in) = []; vert_bias_season(in) = []; vert_bias_refz(in) = []; %remove data from the glacier surface b/c it's a non-stationary feature
            seasonal_bias = [seasonal_bias; vert_bias]; 
            seasonal_id = [seasonal_id; vert_bias_season]; seasonal_refz = [seasonal_refz; vert_bias_refz];
            clear vert_bias* in;
%         end
    end
    writetable(t,csvs(i).name);
    clear t;
end


%%determine whether the seasonal snow- and ice-free residuals to check that
%summer (July-Sept) data should be used, otherwise use spring (April-June) data, using the
%two-sample Kolmogorov-Smirnov test to test for significant differences in distributions
figure; seasonal_colors = [94,60,153; 178,171,210; 253,184,99; 230,97,1]./255; %purple-orange colormap
for j = 1:4
hp(j) = histogram(seasonal_bias(seasonal_id==j & seasonal_refz<=snowline)); hp(j).FaceColor = seasonal_colors(j,:); hold on;
end
xlabel('Elevation difference (m)'); ylabel('Count');
leg = legend(hp,'winter','spring','summer','autumn');
%test difference in spring below snowline & summer residuals
sum_sigtest = kstest2(seasonal_bias(seasonal_id==2 & seasonal_refz<=snowline),seasonal_bias(seasonal_id==3));
if sum_sigtest == 1
    disp('Summer elevation residuals are significantly different than spring residuals below the snowline:');
    fprintf('Summer median +/- MAD errors are %3.2f +/- %3.2f m\n',nanmedian(seasonal_bias(seasonal_id==3)),mad(seasonal_bias(seasonal_id==3),1));
    fprintf('Spring below-snowline median +/- MAD errors are %3.2f +/- %3.2f m\n',...
        nanmedian(seasonal_bias(seasonal_id==2 & seasonal_refz<=snowline)),mad(seasonal_bias(seasonal_id==2 & seasonal_refz<=snowline),1));
end
%user makes the decision to use summer or spring data based on plots & significance testing
prompt = 'Are the summer elevation residuals questionable based on histograms & significance test (y/n)?';
str = input(prompt,'s');
if strmatch(str,'y')==1
    disp('...using spring data for vertical coregistration');
    vertical_adjustment = nanmedian(seasonal_bias(seasonal_id==2 & seasonal_refz<=snowline)); vertical_adjustment(isnan(vertical_adjustment)) = 0;
else
    vertical_adjustment = nanmedian(seasonal_bias(seasonal_id==3)); vertical_adjustment(isnan(vertical_adjustment)) = 0;
end

%vertical coregistration: optional use of transect-specific coregistration
%for ATL08 data if you trust the brightness flag (uncomment optional
%if-else statement within the loop)
disp('Vertically coregistering all transects...');
for i = 1:length(csvs)
    ICESat2 = [csv_path,csvs(i).name]; %compile the file name
    t = readtable(csvs(i).name);
        
    %convert datestamp to decimal dates and add to concatenated file
    yr = str2double(csvs(i).name(7:10));
    if mod(yr,4) == 0
        decidate = (cumdays_leap(str2double(csvs(i).name(11:12)))+str2double(csvs(i).name(13:14)))/sum(modays_leap);
    else
        decidate = (cumdays_norm(str2double(csvs(i).name(11:12)))+str2double(csvs(i).name(13:14)))/sum(modays_norm);
    end
    t.date = repmat(yr+decidate,size(t.Elevation));

    %coregister using median from best snow-free observations identified above
    %if snow is confidently flagged, uncomment if statement for ATL08 to
    %use transect-specific vertical biases during coregistration
%     if contains(csvs(i).name, 'ATL08')
%         if sum(t.Brightness_Flag) > 0
%             t.Elevation_Coregistered = t.Elevation-nanmedian(t.VerticalErrors(t.Brightness_Flag==0));
%             disp([csvs(i).name,' contains ',num2str(sum(t.Brightness_Flag)./height(t)*100),'% snow!']);
%         else %use vertical bias for each profile for snow-free months
%             t.Elevation_Coregistered = t.Elevation-nanmedian(t.VerticalErrors(t.Brightness_Flag==0));
%         end
%     else
        t.Elevation_Coregistered = t.Elevation-vertical_adjustment;
%     end
    
    %compile in one big table
    if i == 1
        writetable(t,[abbrev,'-ICESat2-',acronym,'-params.csv']);
    else
        writetable(t,[abbrev,'-ICESat2-',acronym,'-params.csv'],'WriteMode','Append','WriteVariableNames',false);
    end
    clear t;
end

%calculate elevation bias & RMSE stats after horizontal AND vertical coregistration
disp('Computing elevation residuals (T.differences) from compiled table of all dates');
T = readtable([abbrev,'-ICESat2-',acronym,'-params.csv']);
% [~,T.differences,RMSE] = extract_icesat2_vertical_errors([csv_path,abbrev,'-ICESat2-ATL08-params.csv'],DTM,Ref);
T.differences = T.Elevation_Coregistered - T.ReferenceElevation;
RMSE = sqrt(nansum((T.differences).^2)./length(T.differences));
writetable(T,[abbrev,'-ICESat2-',acronym,'-params.csv']);
%display statistics
disp(['Total regional median +/- MAD elevation difference = ',num2str(nanmedian(T.differences)),...
    '+/-',num2str(mad(T.differences,1)),' (RMSE = ',num2str(RMSE),')']);
disp(['Winter median +/- MAD elevation difference = ',num2str(nanmedian(T.differences(T.season==1))),...
    '+/-',num2str(mad(T.differences(T.season==1),1)),' (RMSE = ',num2str(sqrt(nansum((T.differences(T.season==1)).^2)./length(T.differences(T.season==1)))),')']);
disp(['Spring median +/- MAD elevation difference = ',num2str(nanmedian(T.differences(T.season==2))),...
    '+/-',num2str(mad(T.differences(T.season==2),1)),' (RMSE = ',num2str(sqrt(nansum((T.differences(T.season==2)).^2)./length(T.differences(T.season==2)))),')']);
disp(['Summer median +/- MAD elevation difference = ',num2str(nanmedian(T.differences(T.season==3))),...
    '+/-',num2str(mad(T.differences(T.season==3),1)),' (RMSE = ',num2str(sqrt(nansum((T.differences(T.season==3)).^2)./length(T.differences(T.season==3)))),')']);
disp(['Fall median +/- MAD elevation difference = ',num2str(nanmedian(T.differences(T.season==4))),...
    '+/-',num2str(mad(T.differences(T.season==4),1)),' (RMSE = ',num2str(sqrt(nansum((T.differences(T.season==4)).^2)./length(T.differences(T.season==4)))),')']);
%plot histograms of seasonal differences between ICESat2 & reference elevations
figure;
h(1) = histogram(T.differences); h(1).BinWidth = 1; h(1).FaceColor = [0.5 0.5 0.5]; h(1).EdgeColor = 'none'; hold on;
h(2) = histogram(T.differences(T.season==1)); h(2).BinWidth = 1; h(2).FaceColor = [43,131,186]/255; h(2).FaceAlpha = 1; h(2).EdgeColor = 'k'; h(2).LineWidth = 1.5; hold on;
h(3) = histogram(T.differences(T.season==2)); h(3).BinWidth = 1; h(3).FaceColor = [171,221,164]/255; h(3).FaceAlpha = 0.5; h(3).EdgeColor = 'w'; h(3).LineWidth = 1; hold on;
h(4) = histogram(T.differences(T.season==3)); h(4).BinWidth = 1; h(4).FaceColor = [215,25,28]/255; h(4).FaceAlpha = 0.5; h(4).EdgeColor = 'none'; hold on;
h(5) = histogram(T.differences(T.season==4)); h(5).BinWidth = 1; h(5).FaceColor = [253,174,97]/255; h(5).FaceAlpha = 0.5; h(5).EdgeColor = 'k'; h(5).LineWidth = 1; hold on;
set(gca,'fontsize',16,'xlim',[-15 15]); leg = legend(h,'all','winter','spring','summer','fall');
xlabel('Elevation residuals (m)'); ylabel('Count');
saveas(gcf,[abbrev,'_elevation-residual_histograms.eps'],'epsc');