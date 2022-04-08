# ICESat2-snow-code
Code to analyze ICESat-2 ATL06 and ATL08 data for mountain snow mapping

# Overview
These codes have been designed and tested using ATL08 data for Reynolds Creek Experimental Watershed in Idaho and ATL06 data for Wolverine Glaciers in Alaska. The files should work for any study site for which you have downloaded ATL06 or ATL08 hdf5 files from the National Snow and Ice Data Center (https://nsidc.org/data/icesat-2/products/level-3a). Two site-specific codes are also included that were used for data analysis and figure rendering. These codes can be minimally modified to other sites with the equivalent datasets.

# Installation & Dependencies
These codes were executed using Matlab 2020b and should work with comparable versions. The codes can be downloaded and run for a study site that has (1) high-quality reference elevation map(s) and (2) ICESat-2 level3a hdf5 files. The cmocean color pallettes (https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps) are used for some figures but can be substituted with Matlab's built-in color pallettes. If you find that supporting codes are needed but not provided, please contact the repo owner!

# Execution
Detailed instructions will be provided shortly...
1) write_icesat2_csv: converts downloaded hdf5 files to csvs, extracting basic parameters (x,y,z, etc)

# Outputs
1) Each ICESat-2 hdf5 file (each file is typically one date-stamped track) is converted to a csv file that contain only the most basic variables from the hdf5 file. The x,y,z coordinates in each csv file are adjusted so that the data are coregistered with respect to a reference digital terrain model, with a new csv file produced with the name appended with "-edited". The csv data are organized in a table, with the first row containing descriptive names for the columns.
2) A concatenated csv file is produced for the full dataset for each study site. If terrain parameters have been derived from the reference digital elevation model(s), then terrain parameters will be appended to the study site's csv file. The csv data are organized in a table, with the first row containing descriptive names for the columns.
3) A series of figures can be produced for ATL06 data by modifying the Wolverine_ICESat2_vs_terrain_figures.mlx for a glacierized study site or for ATL08 data by modifying the RCEW_ICESat2_vs_terrain_figures.mlx for a vegetated study site.

# Code Use
This code can be freely used and modified under an MIT license but must be accompanied by a reference to this repo.
