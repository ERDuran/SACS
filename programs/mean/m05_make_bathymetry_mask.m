%% make ETOPO-5 bathymetry mask
% make bathymetry mask
clc
path(pathdef)

% set up main directory
cd /home/z5100174/Desktop/MATLAB

% Add path to the data to be loaded
addpath(genpath('functions'))

% Add path to the data to be loaded
addpath cars_out

% Add path to the gsw TEOS-10 library (including subfolders)
addpath(genpath('teos10_library'))

clear 
% vn is a function that can store the name of a variable
% hence vn(x) returns 'x' as a string.
vn = @(x) inputname(1);

%etopo-5
load raw_data/bathymetry_aus/topog

%
load aus8


%% project topography into 1/8 degree resolution grid
% create four topography masks, one for T and S stations, one for u
% and one for v and one for the divergence (lat_u and lon_v)

% topog_hr_mask = isnan(topog_hr);
% topog_hr(topog_hr_mask) = 0;

% T and S
lat = aus8.temp.var_on_pres.lat;
lon = aus8.temp.var_on_pres.lon;

% we need to smooth out ETOPO-5 because there are some weird discrepancies
smoothing_window = 4;
topog_hr_sm = smooth2a(topog_hr,smoothing_window,smoothing_window);

% T and S interp
topog = round(interp2(...
    lon_topog, lat_topog, topog_hr_sm, lon, lat, 'linear')/10)*10;
topog(isnan(topog)) = 10;

%
topog(topog<-2000) = -2000;


%% plot maps
% p05_plot_CARS_bathymetry


%% save tout ca tout ca
topog_mask.etopo5.topog_orig = topog_hr;
topog_mask.etopo5.topog_smoothed = topog_hr_sm;
topog_mask.etopo5.smoothing_window = smoothing_window;
topog_mask.etopo5.lat = lat_topog;
topog_mask.etopo5.lon = lon_topog;

topog_mask.aus8.topog = topog;
topog_mask.aus8.lat = lat;
topog_mask.aus8.lon = lon;


save('cars_out/topog_mask', 'topog_mask')
disp('topog_mask saved into cars_out !')

