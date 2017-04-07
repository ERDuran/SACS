%% Get aus8 data at selected lon lat window. save into new structure
clc
path(pathdef)

% set up main directory
cd /home/z5100174/Desktop/MATLAB

% Add path to the data to be loaded
addpath cars_grid_mat
addpath(genpath('functions'))

% Add path to the gsw TEOS-10 library (including subfolders)
addpath(genpath('teos10_library'))


%% Load in-situ temperature and practical salinity from aus8
clear

disp('Loading CARS aus8 ...')

data.temp = load('temperature_aus8.mat'); disp('Temperature OK')
data.sal = load('salinity_aus8.mat'); disp('Salinity OK')

data_names = fieldnames(data);

disp('CARS aus8 loaded !')


%% Subset variables of aus8 from lon and lat window
lon = double(data.temp.Var.lon);
lat = double(data.temp.Var.lat);

% Define latitude and longitude ranges to zoom in
lat_range = [-15 -50];
lon_range = [90 155];

% get longitudes within the range
lon_ind = find(lon >= lon_range(1) & lon <= lon_range(2));
% transpose lon to get it as a row vector
lon = lon(lon_ind)';
% get latitudes within the range
lat_ind = find(lat <= lat_range(1) & lat >= lat_range(2));
% flip latitude indexes upside down to make it start from low latitude
lat_ind = flipud(lat_ind);
lat = lat(lat_ind);

for a = 1 : length(data_names)
    % variables
    aus8_raw.(data_names{a}).var.lon = lon;
    aus8_raw.(data_names{a}).var.lat = lat;
    aus8_raw.(data_names{a}).var.depth = double(data.(data_names{a}).Var.depth);
    aus8_raw.(data_names{a}).var.depth_ann = double(data.(data_names{a}).Var.depth_ann);
    aus8_raw.(data_names{a}).var.depth_san = double(data.(data_names{a}).Var.depth_semiann);
    
    % attributes
    aus8_raw.(data_names{a}).att.lon = data.(data_names{a}).var_att.lon;
    aus8_raw.(data_names{a}).att.lat = data.(data_names{a}).var_att.lat;
    aus8_raw.(data_names{a}).att.depth = data.(data_names{a}).var_att.depth;
    aus8_raw.(data_names{a}).att.depth_ann = data.(data_names{a}).var_att.depth_ann;
    aus8_raw.(data_names{a}).att.depth_san = data.(data_names{a}).var_att.depth_semiann;
end


% variables of interest 
Var_names = {...
    'nq', ... % ta ta ta...
    'radius_q', ...
    'mean', ...
    'an_cos', ...
    'an_sin', ...
    'sa_cos', ...
    'sa_sin', ...
    'RMSspatialresid', ...
    'RMSresid', ...
    'sumofwgts'};

for a = 1 : length(data_names)
    for b = 1 : length(Var_names)
        % variable subset
        data_dummy = ...
            data.(data_names{a}).Var.(Var_names{b})(lon_ind,lat_ind,:);
        
        % permute rows and columns
        aus8_raw.(data_names{a}).var.(Var_names{b}) = ...
            permute(data_dummy, [2 1 3]);
        
        % attribute subset
        aus8_raw.(data_names{a}).att.(Var_names{b}) = ...
            data.(data_names{a}).var_att.(Var_names{b});
    end
    
    aus8_raw.(data_names{a}).var_names = Var_names;
end

clear data_dummy lon_ind lat_ind lon_range lat_range
clear lat lon a b
disp('CARS aus8 vars and atts subset and stored into aus8')
save('cars_out/aus8_raw', 'aus8_raw')
disp('CARS aus8_raw subset saved into cars_out')

