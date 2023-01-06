function write_icesat2_csv(inputdir, outputdir, shapefile,starter)
% function WRITE_ICESAT2_CSV reads in an h5 file and outputs csv files of
% individual beams and their pertinent data
% INPUTS: inputdir = directory pointing to .h5 files (ends with '/') (string)
%        outputdir = directory where you want to save csv files (also ends
%                      with '/') (string)
%        shapefile = directory and name of shapefile of the region of
%                      interest that serves for clipping the icesat-2 data
%                      down (string)
%                    Note: the shapefile should be in WGS84 Cartesian Coordinates

% created 21 December 2020 by Colten Elkin (coltenelkin@u.boisestate.edu)
% last modified 22 Feb 2022 by Ellyn Enderlin (ellynenderlin@boisestate.edu)

if(~isfolder(outputdir)) % create the output directory if it doesnt already exist
    mkdir(outputdir)
end

watershed = shaperead(shapefile);
h5files = dir([inputdir,'*.h5']); % pull out the h5 files
beams = {'gt1r', 'gt1l', 'gt2r', 'gt2l', 'gt3r', 'gt3l'}; % list of icesat2 beams for inner loop
for filecounter = starter:length(h5files) % loop through icesat2 files
    % check to see whether it's ATL08 or ATL06
    if strcmp(h5files(filecounter).name(1:5), 'ATL08') == 1 % ATL08 commands
        h5test = [h5files(filecounter).folder,'/',h5files(filecounter).name]; % get string pointing to n-th h5 file
        %create a cell array of dataset names to check that each beam exists
        hinfo = h5info(h5test);
        for j=1:length(hinfo.Groups)
            groups(j) = {hinfo.Groups(j).Name};
        end
        
        %loop through the available beams & pull data
        beamcheck = zeros(size(beams));
        for i = 1:length(groups)
            for j = 1:length(beams)
                if contains(char(groups(i)),char(beams(j)))
                    beamcheck(j) = 1;
                end
            end
        end
        %skip trying to pull data from missing beams
        disp(['Extracting data for ',h5test(end-38:end-3),' (#',num2str(filecounter),')...']);
        for beamcount = 1:length(beams)
            if beamcheck(beamcount) == 1
                beam = beams{beamcount}; disp([beam]); % set beam
                % read in data
                terrain = h5read(h5test,['/',beam,'/land_segments/terrain/h_te_mean']); % read terrain means
                terrain_bestfit = h5read(h5test,['/',beam,'/land_segments/terrain/h_te_best_fit']); % ADDED!!! pull best fit elevations
                lat = h5read(h5test, ['/',beam,'/land_segments/latitude']); % read lats
                lon = h5read(h5test, ['/',beam,'/land_segments/longitude']); % read lons
                bright = h5read(h5test, ['/',beam,'/land_segments/brightness_flag']); % read in brightness of shot
                snow = h5read(h5test, ['/',beam,'/land_segments/segment_snowcover']); % read in if the surface is estimated as 0=water,1=land,2=snow,3=ice
                std = h5read(h5test, ['/',beam,'/land_segments/terrain/h_te_std']); % standard deviation
                can = h5read(h5test, ['/',beam,'/land_segments/canopy/h_canopy']); % canopy elevation
                if ~isempty(terrain)
                    disp('grabbed data');
                end
                
                % crop data to area of interest
                % first, crude cropping
                lonlims = watershed.BoundingBox(:, 1); % save upper and lower longitudes of the watershed
                latlims = watershed.BoundingBox(:, 2); % save upper and lower latitudes of the watershed
                
                % note: this trimming process is ~4x faster than using inpolygon but less precise
                Ix = find(lon > min(lonlims) & lon < max(lonlims)); % find longitudes between limits
                lon = lon(Ix); % cut down based on longitudes
                lat = lat(Ix);
                terrain = terrain(Ix);
                terrain_bestfit = terrain_bestfit(Ix);
                bright = bright(Ix);
                snow = snow(Ix);
                std = std(Ix);
                can = can(Ix);
                
                Iy = find(lat > min(latlims) & lat < max(latlims)); % find lats between limits
                lat = lat(Iy); % cut down based on latitudes
                lon = lon(Iy);
                terrain = terrain(Iy);
                terrain_bestfit = terrain_bestfit(Iy);
                bright = bright(Iy);
                snow = snow(Iy);
                std = std(Iy);
                can = can(Iy);
                
                %final cropping to catchment shapefile bounds
                wats = [watershed.X', watershed.Y']; % save coordinate tuples from the waterhsed shapefile
                Ifinal = inpolygon(lon, lat, wats(:,1), wats(:,2));
                if isempty(Ifinal)
                    disp('...but data are outside ROI');
                end
                lat = lat(Ifinal); % save data as vectors after final clipping
                
                if ~isempty(lat) % continue if data is inside the region of interest
                    lon = lon(Ifinal); % conitnue saving data
                    terrain = terrain(Ifinal);
                    terrain_bestfit = terrain_bestfit(Ifinal);
                    bright = bright(Ifinal);
                    snow = snow(Ifinal);
                    std = std(Ifinal);
                    can = can(Ifinal);
                    can(can > 1000) = nan; % change canopy elevation no data value to nan
                    
                    
                    % use wgs2utm script to write easting and northing values
                    [easts, norths] = wgs2utm(lat, lon);
                    
                    % create a final structure with all of the data
                    s = struct(); % set structure for fields
                    s.Latitude = lat; % set lats
                    s.Longitude = lon; % set lons
                    s.Elevation = terrain; s.Elevation(s.Elevation >= 10e20) = NaN; % set Nans (icesat2 default nan value is 4.028e38)
                    s.Elevation_bestfit = terrain_bestfit; s.Elevation_bestfit(s.Elevation_bestfit >= 10e20) = NaN; % set Nans (icesat2 default nan value is 4.028e38)
                    s.Canopy = can;
                    s.std = std;
                    s.Easting = easts;
                    s.Northing = norths;
                    s.Brightness_Flag = bright;
                    s.Snow_Flag = snow;
                    
                    %write the csv
                    table = struct2table(s); % convert to a table
                    h5filename = h5test(end-38:end-3); % save h5 filename
                    outputname = [outputdir, h5filename, '_', beam, '.csv']; % save full filename
                    writetable(table, outputname)
                    disp(['...output csv for ROI']);
                    clear s table;
                end
                clear terrain lat* lon* I* easts norths bright std snow can terrain_bestfit;
            end
        end
        clear h5test hinfo beamcheck;
        % if not ATL08, ATL06?
    elseif strcmp(h5files(filecounter).name(1:5), 'ATL06') == 1 % enter ATL06 commands
        h5test = [h5files(filecounter).folder,'/',h5files(filecounter).name]; % get string pointing to n-th h5 file
        %create a cell array of dataset names to check that each beam exists
        hinfo = h5info(h5test);
        for j=1:length(hinfo.Groups)
            groups(j) = {hinfo.Groups(j).Name};
        end
        
        %loop through the available beams & pull data
        beamcheck = zeros(size(beams));
        for i = 1:length(groups)
            for j = 1:length(beams)
                if contains(char(groups(i)),char(beams(j)))
                    beamcheck(j) = 1;
                end
            end
        end
        %skip trying to pull data from missing beams
        disp(['Extracting data for ',h5test(end-38:end-3),' (#',num2str(filecounter),')...']);
        for beamcount = 1:length(beams)
            if beamcheck(beamcount) == 1
                beam = beams{beamcount}; disp([beam]); % set beam
                % read in data
                terrain = h5read(h5test, ['/',beam,'/land_ice_segments/h_li']); % read terrain means
                lat = h5read(h5test,  ['/',beam,'/land_ice_segments/latitude']); % read lats
                lon = h5read(h5test,  ['/',beam,'/land_ice_segments/longitude']); % read lons
                geo_err = h5read(h5test,  ['/',beam,'/land_ice_segments/sigma_geo_h']); % read in vertical geolocation error
                std = h5read(h5test,  ['/',beam,'/land_ice_segments/h_li_sigma']); % standard deviation
                geoid = h5read(h5test,  ['/',beam,'/land_ice_segments/dem/geoid_h']); % standard deviation
                if ~isempty(terrain)
                    disp('grabbed data');
                end
                
                % crop data to area of interest
                % first, crude cropping
                lonlims = watershed.BoundingBox(:, 1); % save upper and lower longitudes of the watershed
                latlims = watershed.BoundingBox(:, 2); % save upper and lower latitudes of the watershed
                
                
                % note: this trimming process is ~4x faster than using inpolygon but less precise
                Ix = find(lon > min(lonlims) & lon < max(lonlims)); % find longitudes between limits
                lon = lon(Ix); % cut down based on longitudes
                lat = lat(Ix);
                terrain = terrain(Ix);
                geo_err = geo_err(Ix);
                std = std(Ix);
                geoid = geoid(Ix);
                
                Iy = find(lat > min(latlims) & lat < max(latlims)); % find lats between limits
                lat = lat(Iy); % cut down based on latitudes
                lon = lon(Iy);
                terrain = terrain(Iy);
                geo_err = geo_err(Iy);
                std = std(Iy);
                geoid = geoid(Iy);
                
                %final cropping to catchment shapefile bounds
                wats = [watershed.X', watershed.Y']; % save coordinate tuples from the waterhsed shapefile
                Ifinal = inpolygon(lon, lat, wats(:,1), wats(:,2));
                if isempty(Ifinal)
                    disp('...but data are outside ROI');
                end
                lat = lat(Ifinal); % save data as vectors after final clipping
                
                if ~isempty(lat) % continue if data is inside the region of interest
                    lon = lon(Ifinal); % conitnue saving data
                    terrain = terrain(Ifinal);
                    geo_err = geo_err(Ifinal);
                    std = std(Ifinal);
                    geoid = geoid(Ifinal);
                    
                    % use wgs2utm script to write easting and northing values
                    [easts, norths] = wgs2utm(lat, lon);
                    
                    % create a final structure with all of the data
                    s = struct(); % set structure for fields
                    s.Latitude = lat; % set lats
                    s.Longitude = lon; % set lons
                    s.Elevation = terrain; s.Elevation(s.Elevation >= 10e20) = NaN; % set Nans (icesat2 default nan value is 4.028e38)
                    s.std = std;
                    s.Easting = easts;
                    s.Northing = norths;
                    s.Vert_Geo_error = geo_err;
                    s.Geoid = geoid;
                    
                    %write the csv
                    table = struct2table(s); % convert to a table
                    h5filename = h5test(end-38:end-3); % save h5 filename
                    outputname = [outputdir, h5filename, '_', beam, '.csv']; % save full filename
                    writetable(table, outputname)
                    disp(['...output csv for ROI']);
                    clear s table;
                end
                clear terrain lat* lon* I* easts norths geo_err std geoid;
            end
        end
        clear h5test hinfo beamcheck;
    end
end
