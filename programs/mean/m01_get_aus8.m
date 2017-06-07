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
aus8.lat = lat(lat_ind);

% get longitudes within the range
lon_ind = find(lon >= lon_range(1) & lon <= lon_range(2));
% transpose lon to get it as a row vector
aus8.lon = lon(lon_ind)';

depth_ind = find(depth >= depth_max);


aus8.depth = depth(depth_ind);
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


%%
save([data_path 'SACS_data/aus8'], 'aus8')
disp('aus8 DONE')

