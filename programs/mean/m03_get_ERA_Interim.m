%% get wind data
clc
path(pathdef)

% set up main directory
cd /home/z5100174/Desktop/MATLAB

% Add path of useful functions (including subfolders)
addpath(genpath('functions'))

% Add path of the calculations output folder
addpath raw_data/wind

clear

[var, gen_att, var_att] = ...
    nc2mat('_grib2netcdf-atls13-a562cefde8a29a7288fa0b8b7f9413f7-9jzuJX.nc', ...
    'ALL');


%% get vars
lat = double(var.latitude);
lon = double(var.longitude)';
%
time_raw = double(var.time);
time = datestr(datenum(1900,1,1) + time_raw/24);

%
inss_raw = var.inss;
inss = permute(inss_raw, [2, 1, 3]);

inss_annual = mean(inss,3);

%
iews_raw = var.iews;
iews = permute(iews_raw, [2, 1, 3]);

iews_annual = mean(iews,3);

%
ERA_Interim.lat = lat;
ERA_Interim.lon = lon;

ERA_Interim.time.units = var_att.time.units;
ERA_Interim.time.time_num = time_raw;
ERA_Interim.time.time_str = time;

ERA_Interim.inss.units = var_att.inss.units;
ERA_Interim.inss.long_name = var_att.inss.long_name;
ERA_Interim.inss.value = inss;
ERA_Interim.inss.annual = inss_annual;

ERA_Interim.iews.units = var_att.iews.units;
ERA_Interim.iews.long_name = var_att.iews.long_name;
ERA_Interim.iews.value = iews;
ERA_Interim.iews.annual = iews_annual;


%% save
save('cars_out/ERA_Interim', 'ERA_Interim')
disp('ERA_Interim saved in cars_out/ERA_Interim !')


