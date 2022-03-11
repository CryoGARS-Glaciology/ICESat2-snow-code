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
else
    disp('If using DEM time series, make sure they are vertically coregistered (done in concatenate_glacier_DEM_timeseries.m)');
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

%% Coregistration
disp('Horizontal coregistration. This takes a LONG LONG time!');

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

%look for the coregistration output files & reload them if they exist
cd_to_csv = ['cd ',csv_path]; eval(cd_to_csv);
output_shifts = dir([abbrev,'-Abest.mat']);
if ~isempty(output_shifts)
   load([abbrev,'-Abest.mat']); load([abbrev,'-RMSDbest.mat']); 
end

%transect csvs
cd_to_csv = ['cd ',csv_path]; eval(cd_to_csv);
csvs = dir([acronym,'*.csv']); %if running more than once, rename original csvs to '*raw.csv' and change the search here to find files with that ending
for i = 1:length(csvs)
    filenamelength(i) = length(csvs(i).name);
end
csvs(filenamelength>45) = []; %automatically ignore edited files (longer names)
clear filenamelength;

%eliminate all short individual transects
for i = 1:length(csvs)
    T = readtable(csvs(i).name);
    
    %check number of observations & delete all files from that date if
    %there are <3 observations (minimum needed to determine footprint)
    if length(T.Elevation) < 3
        delete(csvs(i).name);
    end
end

%identify the raw ICESat-2 data csv files (looped to filter bad dates)
for k = 1:2
    cd_to_csv = ['cd ',csv_path]; eval(cd_to_csv);
    csvs = dir([acronym,'*.csv']); %if running more than once, rename original csvs to '*raw.csv' and change the search here to find files with that ending
    for i = 1:length(csvs)
        filenamelength(i) = length(csvs(i).name);
    end
    csvs(filenamelength>45) = []; %automatically ignore edited files (longer names)
    clear filenamelength;
    
    %identify unique dates from csv filesnames (each date can have up to 6 csvs)
    for i = 1:length(csvs)
        csv_dates(i,:) = csvs(i).name(7:20);
    end
    [unique_dates,unique_refs,unique_inds] = unique(csv_dates,'rows');
    clear csv_dates;
    
    %check that you have enough good data for each data (delete if <3 points)
    if k == 1
        %concatenated dates
        for i = 1:length(unique_refs)
            %find all files with the same date
            filerefs = find(unique_inds==i);
            
            %concatenate all files with the same date
            t = [];
            for j = 1:length(filerefs)
                T = readtable(csvs(filerefs(j)).name); t = [t; T]; clear T;
            end
            
            %check number of observations & delete all files from that date if
            %there are <3 observations (minimum needed to determine footprint)
            if length(t.Elevation) <3
                for j = 1:length(filerefs)
                    delete(csvs(filerefs(j)).name);
                end
            end
        end
    end
end

%loop through the csvs & determine the horizontal offsets
for i = 1:length(unique_refs)
    tic;
    disp(['Coregistering data from date #',num2str(i),' of ',num2str(length(unique_refs)),' (',unique_dates(i,1:4),'/',unique_dates(i,5:6),'/',unique_dates(i,7:8),')']);
    %identify
    filerefs = find(unique_inds==i);
    for j = 1:length(filerefs)
        ICESat2_files(j,:) = [csvs(filerefs(j)).folder,'/',csvs(filerefs(j)).name]; %compile the file name
    end
    
    %if analyzing ATL06 data for a glacier, select the most recent DEM
    if contains(acronym,'ATL06')

        %convert datestamp to decimal dates
        yr = str2double(csvs(filerefs(1)).name(7:10));
        if mod(yr,4) == 0
            decidate = (cumdays_leap(str2double(csvs(filerefs(1)).name(11:12)))+str2double(csvs(filerefs(1)).name(13:14)))/sum(modays_leap);
        else
            decidate = (cumdays_norm(str2double(csvs(filerefs(1)).name(11:12)))+str2double(csvs(filerefs(1)).name(13:14)))/sum(modays_norm);
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
    myFuncHandle = @(A)coregister_icesat2(ICESat2_files,DTM,Ref,A); %create the handle to call the coregistration function
    [Abest(i,:),RMSDbest(i)] = fminsearch(myFuncHandle,[0,0]); %initial horizontal offset estimate = [0,0] = [0 m East, 0 m North]
    save([abbrev,'-Abest.mat'],'Abest','-v7.3'); save([abbrev,'-RMSDbest.mat'],'RMSDbest','-v7.3');
    fprintf('x-offset = %5.2f m & y-offset = %5.2f m w/ RMSD = %5.2f m \n',Abest(i,1),Abest(i,2),RMSDbest(i));
    clear vals;
    toc; %display processing time
    
    %horizontally coregister & write the reference DTM and vertical offset values to the csv file
    for j = 1:length(filerefs)
        %read the data
        T = readtable(csvs(filerefs(j)).name);
        %specify variables with column indices in square brackets if you want to cull data
        if contains(csvs(filerefs(j)).name, 'ATL08')
            t = T(:,:);
        else
            t = T(:,:);
        end
        %apply horizontal coregistration shift that minimizes the rmsd in elevations determined using fminsearch above
        if ~isnan(Abest(i,1))
            t.Northing = t.Northing+Abest(i,2); t.Easting = t.Easting+Abest(i,1);
        end
        writetable(t,[csvs(filerefs(j)).name(1:end-4),'-edited.csv']);
        
        %extract reference elevations and vertical errors
        if length(t.Elevation) > 1
            [t.ReferenceElevation,t.VerticalErrors,~] = extract_icesat2_vertical_errors([csvs(filerefs(j)).name(1:end-4),'-edited.csv'],DTM,Ref);
            fprintf('Median vertical bias after horizontal coregistration = %5.2f m \n',nanmedian(t.VerticalErrors));
        else
            t.ReferenceElevation = NaN; t.VerticalErrors = NaN;
            disp('Too few data to coregister!');
        end
        writetable(t,[csvs(filerefs(j)).name(1:end-4),'-edited.csv']);
        
        clear T t;
    end
    disp(['Adjusted coordinates and reference elevations saved for ',unique_dates(i,1:4),'/',unique_dates(i,5:6),'/',unique_dates(i,7:8)]);
    disp(' '); %leave a space in the command line after the outputs
end

%% apply median snow-free offset to vertically coregister 
close all; disp('Vertical coregistration:');

%read in the snow-off reference elevation map
cd_to_DTM = ['cd ',DTM_path]; eval(cd_to_DTM);
if contains(DTM_name,'.tif')
    [DTM,Ref] = readgeoraster(DTM_name);
    yr = str2double(DTM_date(1:4));
    if mod(yr,4) == 0
        decidate = (cumdays_leap(str2double(DTM_date(5:6)))+str2double(DTM_date(7:8)))/sum(modays_leap);
    else
        decidate = (cumdays_norm(str2double(DTM_date(5:6)))+str2double(DTM_date(7:8)))/sum(modays_norm);
    end
    deciyear = yr+decidate;
    DEMdate = deciyear; clear yr decidate deciyear;
elseif contains(DTM_name,'.mat')
    load_DTM = ['load ',DTM_name]; eval(load_DTM);
    for i = 1:length(Z)
        DEMdate(i) = Z(i).deciyear;
    end
end

%identify the ICESat-2 data csv files
cd_to_csv = ['cd ',csv_path]; eval(cd_to_csv);
csvs = dir('*-edited.csv');

%loop through the csvs, assign seasons, & compile into one table
seasonal_zbias = []; seasonal_id = []; seasonal_yrref = []; seasonal_zref = [];
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

	%if using multiple time-stamped DEMs, identify the DEM used
    if length(DEMdate) > 1
        %convert csv datestamp to decimal dates
        yr = str2double(csvs(i).name(7:10));
        if mod(yr,4) == 0
            decidate = (cumdays_leap(str2double(csvs(i).name(11:12)))+str2double(csvs(i).name(13:14)))/sum(modays_leap);
        else
            decidate = (cumdays_norm(str2double(csvs(i).name(11:12)))+str2double(csvs(i).name(13:14)))/sum(modays_norm);
        end
        deciyear = yr+decidate;
        
        %identify the most recent DEM
        yr_ref = find(DEMdate<=deciyear,1,'last');
        t.DEMdate = repmat(DEMdate(yr_ref),size(t.season));
        clear yr_ref;
        clear yr decidate deciyear;
    else
        t.DEMdate = repmat(DEMdate,size(t.season));
    end
    
    %compile seasonal elevation biases
    if contains(csvs(i).name, 'ATL08')
        seasonal_zbias = [seasonal_zbias; t.VerticalErrors(t.Brightness_Flag==0)]; %make sure only snow-free elevation differences are used
        seasonal_id = [seasonal_id; t.season(t.Brightness_Flag==0)];
        seasonal_yrref = [seasonal_yrref; t.DEMdate(t.Brightness_Flag==0)];
        seasonal_zref = [seasonal_zref; t.ReferenceElevation(t.Brightness_Flag==0)];
    else
        vert_bias = t.VerticalErrors(t.ReferenceElevation<=snowline); %must exclude data above snowline year-round for the glacier since the snowline marks the elevation of year-round snow
        vert_bias_season = t.season(t.ReferenceElevation<=snowline); vert_bias_yr = t.DEMdate(t.ReferenceElevation<=snowline);
        vert_bias_refz = t.ReferenceElevation(t.ReferenceElevation<=snowline);
        in = inpolygon(t.Easting(t.ReferenceElevation<=snowline), t.Northing(t.ReferenceElevation<=snowline), S.X, S.Y); %identify footprints in the glacier outline
        vert_bias(in) = []; vert_bias_season(in) = []; vert_bias_yr(in) = []; vert_bias_refz(in) = []; %remove data from the glacier surface b/c it's a non-stationary feature
        seasonal_zbias = [seasonal_zbias; vert_bias];
        seasonal_id = [seasonal_id; vert_bias_season]; 
        seasonal_yrref = [seasonal_yrref; vert_bias_yr];
        seasonal_zref = [seasonal_zref; vert_bias_refz];
        clear vert_bias* in;
    end
    writetable(t,csvs(i).name);
    clear t;
end


%%determine whether the seasonal snow- and ice-free residuals to check that
%summer (July-Sept) data should be used, otherwise use spring (April-June) data, using the
%two-sample Kolmogorov-Smirnov test to test for significant differences in distributions
figure; set(gcf,'position',[50 50 800 1000]); clear hp;
seasonal_colors = [43,131,186; 171,221,164; 215,25,28; 253,174,97]./255;
for i = 1:length(DEMdate)
    subplot(length(DEMdate),1,i); subpos = get(gca,'position');
    for j = [4 1 2 3];
        hp(j) = histogram(seasonal_zbias(seasonal_id==j & seasonal_zref<=snowline & seasonal_yrref == DEMdate(i)),[floor(nanmedian(seasonal_zbias)-3*1.4826*mad(seasonal_zbias,1)):1:ceil(nanmedian(seasonal_zbias)+3*1.4826*mad(seasonal_zbias,1))],'Normalization','pdf'); hp(j).FaceColor = seasonal_colors(j,:); hold on;
    end
    set(gca,'fontsize',16); ylims = get(gca,'ylim');
    if length(DEMdate) > 1
        text(floor(nanmedian(seasonal_zbias)-3*1.4826*mad(seasonal_zbias,1)),max(ylims)-0.05*range(ylims),[' Year: ',num2str(floor(DEMdate(i))),'-',num2str(ceil(DEMdate(i)))],'fontsize',14);
    end
if i == 1;
title('Off-ice elevation differences before vertical coregistration');
end
end
xlabel('Elevation difference (m)','fontsize',16); ylabel('Probability density','fontsize',16);
leg = legend(hp,'winter','spring','summer','autumn');
%test difference in spring below snowline & summer residuals
sum_sigtest = kstest2(seasonal_zbias(seasonal_id==2 & seasonal_zref<=snowline),seasonal_zbias(seasonal_id==3));
if sum_sigtest == 1
    disp('Summer elevation residuals are significantly different than spring residuals below the snowline:');
    fprintf('Summer median +/- MAD errors are %3.2f +/- %3.2f m\n',nanmedian(seasonal_zbias(seasonal_id==3)),mad(seasonal_zbias(seasonal_id==3),1));
    fprintf('Spring below-snowline median +/- MAD errors are %3.2f +/- %3.2f m\n',...
        nanmedian(seasonal_zbias(seasonal_id==2 & seasonal_zref<=snowline)),mad(seasonal_zbias(seasonal_id==2 & seasonal_zref<=snowline),1));
end
%user makes the decision to use summer or spring data based on plots & significance testing
prompt = 'Are the summer elevation residuals questionable based on histograms & significance test (y/n)?';
str = input(prompt,'s');
if strmatch(str,'y')==1
    disp('...using spring data for vertical coregistration');
    vertical_adjustment = nanmedian(seasonal_zbias(seasonal_id==2 & seasonal_zref<=snowline)); vertical_adjustment(isnan(vertical_adjustment)) = 0;
else
    disp('...using summer data for vertical coregistration');
    vertical_adjustment = nanmedian(seasonal_zbias(seasonal_id==3)); vertical_adjustment(isnan(vertical_adjustment)) = 0;
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
    
    
    %extract beam info from the file name and add to concatenated file
    t.beam = repmat(csvs(i).name(end-14:end-11),size(t.Elevation));
    

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
figure; clear h;
h(1) = histogram(T.differences(T.season==3),'Normalization','pdf'); h(1).BinWidth = 0.5; h(1).FaceColor = seasonal_colors(3,:); h(1).FaceAlpha = 1; h(1).EdgeColor = seasonal_colors(3,:); hold on;
h(2) = histogram(T.differences(T.season==4),'Normalization','pdf'); h(2).BinWidth = 0.5; h(2).FaceColor = seasonal_colors(4,:); h(2).FaceAlpha = 0.5; h(2).EdgeColor = seasonal_colors(4,:); h(2).LineWidth = 1; hold on;
h(3) = histogram(T.differences(T.season==1),'Normalization','pdf'); h(3).BinWidth = 0.5; h(3).FaceColor = seasonal_colors(1,:); h(3).FaceAlpha = 0.5; h(3).EdgeColor = seasonal_colors(1,:); h(3). LineStyle = '--'; hold on;
h(4) = histogram(T.differences(T.season==2),'Normalization','pdf'); h(4).BinWidth = 0.5; h(4).FaceColor = seasonal_colors(2,:); h(4).FaceAlpha = 0.5; h(4).EdgeColor = seasonal_colors(2,:); h(4).LineWidth = 1; h(4). LineStyle = ':'; hold on;
set(gca,'fontsize',16,'xlim',[floor(nanmedian(T.differences)-3*1.4826*mad(T.differences,1)) ceil(nanmedian(T.differences)+3*1.4826*mad(T.differences,1))]); 
leg = legend(h,'summer','autumn','winter','spring');
xlabel('Elevation residuals (m)'); ylabel('Probability density');
saveas(gcf,[abbrev,'_elevation-residual_histograms.eps'],'epsc');
