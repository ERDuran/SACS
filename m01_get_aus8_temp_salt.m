%% Get aus8 data at selected lon lat depth
clearvars('-except', '*_path')

% temperature
[temp_data, temp_gen_att, temp_att] = ...
    nc2mat([data_path ...
    'CARS/aus8/temperature_aus8_2013.nc'], ...
    'ALL');

% salinity
[salt_data, salt_gen_att, salt_att] = ...
    nc2mat([data_path ...
    'CARS/aus8/salinity_aus8_2013.nc'], ...
    'ALL');

lat = double(temp_data.lat);
lon = double(temp_data.lon);
depth = -double(temp_data.depth);
depth_ann = -double(temp_data.depth_ann);
depth_semiann = -double(temp_data.depth_semiann);

temp.mean = temp_data.mean;
temp.nq = temp_data.nq;
temp.radius_q = temp_data.radius_q;
temp.an_cos = temp_data.an_cos;
temp.an_sin = temp_data.an_sin;
temp.sa_cos = temp_data.sa_cos;
temp.sa_sin = temp_data.sa_sin;

salt.mean = salt_data.mean;
salt.nq = salt_data.nq;
salt.radius_q = salt_data.radius_q;
salt.an_cos = salt_data.an_cos;
salt.an_sin = salt_data.an_sin;
salt.sa_cos = salt_data.sa_cos;
salt.sa_sin = salt_data.sa_sin;


%% Subset variables of aus8 from lon and lat window
% Define latitude and longitude ranges to zoom in
lat_range = [-30 -50];
lon_range = [108 153];
depth_max = -2000;

% get latitudes within the range
lat_ind = find(lat <= lat_range(1) & lat >= lat_range(2));
% flip latitude indexes upside down
lat_ind = flipud(lat_ind);

% get longitudes within the range
lon_ind = find(lon >= lon_range(1) & lon <= lon_range(2));

depth_ind = find(depth >= depth_max);

aus8_coor.lat = lat(lat_ind);
aus8_coor.lon = lon(lon_ind)';
aus8_coor.depth = depth(depth_ind);
save([data_path 'SACS_data/aus8_coor'], 'aus8_coor')

aus8.depth_ann = depth_ann;
aus8.depth_semiann = depth_semiann;

aus8.temp.mean = permute(temp.mean(lon_ind,lat_ind,depth_ind),[2,1,3]);
aus8.temp.nq = permute(temp.nq(lon_ind,lat_ind,depth_ind),[2,1,3]);
aus8.temp.radius_q = ...
    permute(temp.radius_q(lon_ind,lat_ind,depth_ind),[2,1,3]);
aus8.temp.an_cos = permute(temp.an_cos(lon_ind,lat_ind,:),[2,1,3]);
aus8.temp.an_sin = permute(temp.an_sin(lon_ind,lat_ind,:),[2,1,3]);
aus8.temp.sa_cos = permute(temp.sa_cos(lon_ind,lat_ind,:),[2,1,3]);
aus8.temp.sa_sin = permute(temp.sa_sin(lon_ind,lat_ind,:),[2,1,3]);

aus8.salt.mean = permute(salt.mean(lon_ind,lat_ind,depth_ind),[2,1,3]);
aus8.salt.nq = permute(salt.nq(lon_ind,lat_ind,depth_ind),[2,1,3]);
aus8.salt.radius_q = ...
    permute(salt.radius_q(lon_ind,lat_ind,depth_ind),[2,1,3]);
aus8.salt.an_cos = permute(salt.an_cos(lon_ind,lat_ind,:),[2,1,3]);
aus8.salt.an_sin = permute(salt.an_sin(lon_ind,lat_ind,:),[2,1,3]);
aus8.salt.sa_cos = permute(salt.sa_cos(lon_ind,lat_ind,:),[2,1,3]);
aus8.salt.sa_sin = permute(salt.sa_sin(lon_ind,lat_ind,:),[2,1,3]);


% make sure temp and salt have the same dry points
temp_nan = isnan(aus8.temp.mean);
salt_nan = isnan(aus8.salt.mean);
aus8.temp.mean(salt_nan) = NaN;
aus8.salt.mean(temp_nan) = NaN;


%% pcolor_maker 
fig_n = 1;
rowcols = [3 2];
rowcols_size = [8 5]; % cm
margs = [1 1 1 1]; % cm
gaps = [1 1]; % cm

for sp = 1 : rowcols(1)*rowcols(2)
    axis_setup{sp} = ...
        [aus8_coor.lon(1) aus8_coor.lon(end) ...
        aus8_coor.lat(end) aus8_coor.lat(1)];
    x{sp} = aus8_coor.lon;
    y{sp} = aus8_coor.lat;
end
    
z_ind = {0 0 -400 -400 -1000 -1000};
data = {...
    aus8.temp.mean(:,:,depth==z_ind{1}), ...
    aus8.salt.mean(:,:,depth==z_ind{2}), ...
    aus8.temp.mean(:,:,depth==z_ind{3}), ...
    aus8.salt.mean(:,:,depth==z_ind{4}), ...
    aus8.temp.mean(:,:,depth==z_ind{5}), ...
    aus8.salt.mean(:,:,depth==z_ind{6})};
minmax = {...
    [8 20], ...
    [34 36], ...
    [min(min(data{3}))+2 max(max(data{3}))-2], ...
    [min(min(data{4}))+0.2 max(max(data{4}))-0.2], ...
    [min(min(data{5})) max(max(data{5}))], ...
    [min(min(data{6})) max(max(data{6}))]};
cmaps_levels = 12;
cmaps = {...
    flipud(othercolor('RdYlBu8', cmaps_levels)), ...
    flipud(othercolor('PuOr8', cmaps_levels)), ...
    flipud(othercolor('RdYlBu8', cmaps_levels)), ...
    flipud(othercolor('PuOr8', cmaps_levels)), ...
    flipud(othercolor('RdYlBu8', cmaps_levels)), ...
    flipud(othercolor('PuOr8', cmaps_levels))};
titles = {...
    ['aus8 mean $temp$ $z=' num2str(z_ind{1}) '$ ($^{\circ}C$)'], ...
    ['aus8 mean $salt$ $z=' num2str(z_ind{2}) '$ ($psu$)'], ...
    ['aus8 mean $temp$ $z=' num2str(z_ind{3}) '$ ($^{\circ}C$)'], ...
    ['aus8 mean $salt$ $z=' num2str(z_ind{4}) '$ ($psu$)'], ...
    ['aus8 mean $temp$ $z=' num2str(z_ind{5}) '$ ($^{\circ}C$)'], ...
    ['aus8 mean $salt$ $z=' num2str(z_ind{6}) '$ ($psu$)']};
font_size = 9;
fig_color = [0.7 0.7 0.7];

fig = pcolor_maker(...
    fig_n, rowcols, rowcols_size, margs, gaps, ...
    x, y, data, axis_setup, minmax, cmaps, ...
    titles, font_size, fig_color);
outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
export_fig(fig, ...
    [figures_path mfilename '/' scriptname(1:3) ...
    '_fig' num2str(fig_n) '_'], ...
    '-m3')
close


%%
save([data_path 'SACS_data/aus8'], 'aus8')
disp('aus8 DONE')

