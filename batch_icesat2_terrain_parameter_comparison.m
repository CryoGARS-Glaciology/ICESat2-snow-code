%%%
%Loop through all compiled csv file containing ICESat-2 ATL08 or ATL06 data
%and extract terrain parameters.
%1) Path where you keep your codes 
%2) Path for the DTM
%3) Name of the DTM
%4) Path for the csv file (set-up for the table of compiled data for all dates)
%5) Site abbreviation used for compiled outputs
%6) ICESat-2 data acronym (e.g., ATL06, ATL08)
%7) Path with filename for glacier outline (if working with ATL06 data)
%8) Parameter files (set-up for slope, aspect, & ruggedness geotiffs
%produced from the reference DTM using QGIS build-in functions)

%Outputs: 
%1) Appended csv file that includes the specified terrain parameter for
%each footprint


%%%
%% Inputs
clearvars; close all; 

%path for the code
addpath('/Users/ellynenderlin/Research/NASA_CryoIdaho/ICESat2-snow-code/')

%terrain parameter file (be sure the path ends in a /)
TP_path = '/Users/ellynenderlin/Research/NASA_CryoIdaho/mountains/RCEW/DEMs/';
TP_name = 'RCEW_1m_WGS84UTM11_WGS84.tif'; %list the DEM name and it will find all associated terrain param files

%site concanated csv (be sure the path ends in a /)
csv_path = '/Users/ellynenderlin/Research/NASA_CryoIdaho/mountains/RCEW/csvs/';

%site abbreviation for file names
abbrev = 'RCEW'; 

%ICESat-2 product acronym
acronym = 'ATL08';

%recommended terrain parameters to add to the table
if contains(acronym,'ATL08')
    terrain_params = {'slope','aspect','ruggedness','VegHeight'};
else
    terrain_params = {'slope','aspect','ruggedness'};
end

% %if ATL06: glacier outline
% if contains(acronym,'ATL06')
%     glims = shaperead('/Users/ellynenderlin/Research/NASA_CryoIdaho/Colten-prelim-data/glims_polygons_2d.shp'); % load in polygon
%     [shpx,shpy,utmzone] = wgs2utm(glims.Y,glims.X);
% end


%% 
disp('This may take a while, especially if using time-stamped reference data!');

%load the ICESat-2 data table
cd_to_csv = ['cd ',csv_path]; eval(cd_to_csv);
T = readtable([abbrev,'-ICESat2-',acronym,'-params.csv']);
icesat2 = [csv_path,abbrev,'-ICESat2-',acronym,'-params.csv'];

%determine whether this is a data update or a new data grab
prompt = 'Are you replacing terrain parameters with updated data (y/n)?';
str = input(prompt,'s');

%pull metrics for each terrain parameter
for i = 1:length(terrain_params)
    
    %load the terrain parameters & extract metrics
    cd_to_DTM = ['cd ',TP_path]; eval(cd_to_DTM);
    
    %if adding the terrain parameter for the first time
    if contains(str,'n')
        if ~ismember(char(terrain_params(i)), T.Properties.VariableNames)
            if contains(TP_name,'.tif')
                tifs = dir('*.tif');
                for j = 1:length(tifs)
                    if contains(tifs(j).name,char(terrain_params(i)))
                        disp(['Extracting stats for ',char(terrain_params(i))]);
                        [tif,R2] = readgeoraster(tifs(j).name);
                        [params,~,~] = calc_icesat2_params(icesat2, tif, R2);
                        params(params==-9999) = NaN;
                        T = addvars(T,params,'NewVariableNames',char(terrain_params(i)));
                        clear params tif R2;
                    end
                end
            elseif contains(TP_name,'.mat')
                mats = dir('*.mat');
                for j = 1:length(mats)
                    if contains(mats(j).name,char(terrain_params(i)))
                        disp(['Extracting stats for ',char(terrain_params(i))]);
                        load(mats(j).name); %load the concatenated matfile
                        [params,~,~] = calc_icesat2_params(icesat2, Z);
                        params(params==-9999) = NaN;
                        T = addvars(T,params,'NewVariableNames',char(terrain_params(i)));
                        clear params tif R2;
                    end
                end
            end
        end
    else
        if contains(TP_name,'.tif')
                tifs = dir('*.tif');
                for j = 1:length(tifs)
                    if contains(tifs(j).name,char(terrain_params(i)))
                        disp(['Extracting stats for ',char(terrain_params(i))]);
                        [tif,R2] = readgeoraster(tifs(j).name);
                        [params,~,~] = calc_icesat2_params(icesat2, tif, R2);
                        params(params==-9999) = NaN;
                        replace_data = ['T.',char(terrain_params(i)),' = params;']; eval(replace_data);
                        clear params tif R2;
                    end
                end
            elseif contains(TP_name,'.mat')
                mats = dir('*.mat');
                for j = 1:length(mats)
                    if contains(mats(j).name,char(terrain_params(i)))
                        disp(['Extracting stats for ',char(terrain_params(i))]);
                        load(mats(j).name); %load the concatenated matfile
                        [params,~,~] = calc_icesat2_params(icesat2, Z);
                        params(params==-9999) = NaN;
                        replace_data = ['T.',char(terrain_params(i)),' = params;']; eval(replace_data);
                        clear params tif R2;
                    end
                end
            end
    end
    cd_to_csv = ['cd ',csv_path]; eval(cd_to_csv);
    writetable(T,[abbrev,'-ICESat2-',acronym,'-params.csv']);
    disp('Resaved table');
end
disp('Done extracting terrain parameters!');
