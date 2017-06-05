%% make ETOPO-5 bathymetry mask
clearvars('-except', '*_path')

%
[topo_data, topo_gen_att, topo_att] = ...
    nc2mat([data_path ...
    'ETOPO5/etopo5.nc'], ...
    'ALL');

lat = flipud(topo_data.topo_lat);
lon = topo_data.topo_lon';
topo = flipud(double(topo_data.topo)');


%%
lat_range = [-30 -50];
lon_range = [108 153];

lat_aus8 = (lat_range(1) : -1/8 : lat_range(2))';
lon_aus8 = lon_range(1) : 1/8 : lon_range(2);

% get latitudes within the range
lat_ind = find(lat <= lat_range(1) & lat >= lat_range(2));
% flip latitude indexes upside down
ETOPO5.lat = lat(lat_ind);

% get longitudes within the range
lon_ind = find(lon >= lon_range(1) & lon <= lon_range(2));
% transpose lon to get it as a row vector
ETOPO5.lon = lon(lon_ind)';

topo_orig = topo(lat_ind, lon_ind);
topo_orig_interp = interp2(...
    ETOPO5.lon, ETOPO5.lat, topo_orig, ...
    lon_aus8, lat_aus8, 'linear');

% we need to smooth out ETOPO-5 because there are some weird discrepancies
ETOPO5.sm_factor = 1;
ETOPO5.sm_window = (2*ETOPO5.sm_factor+1)*1/12;
topo_sm = smooth2a(topo_orig,ETOPO5.sm_factor);

% T and S interp
topo_sm_interp = interp2(...
    ETOPO5.lon, ETOPO5.lat, topo_sm, lon_aus8, lat_aus8, 'linear');

%
topo_orig(topo_orig>0) = 0;
topo_orig_interp(topo_orig_interp>0) = 0;
topo_sm(topo_sm>0) = 0;
topo_sm_interp(topo_sm_interp>0) = 0;

%
topo_orig_interp = round(topo_orig_interp/10)*10;
topo_sm_interp = round(topo_sm_interp/10)*10;


%%
ETOPO5.topo_orig = topo_orig;
ETOPO5.topo_orig_interp = topo_orig_interp;
ETOPO5.topo_sm = topo_sm;
ETOPO5.topo_sm_interp = topo_sm_interp;


%%
save([data_path 'SACS_data/ETOPO5'], 'ETOPO5')
disp('ETOPO5 DONE')

