%% make ETOPO-5 bathymetry mask
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
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
% topo_orig_interp = round(topo_orig_interp/10)*10;
% topo_sm_interp = round(topo_sm_interp/10)*10;

%
ETOPO5.topo_orig = topo_orig;
ETOPO5.topo_orig_interp = topo_orig_interp;
ETOPO5.topo_sm = topo_sm;
ETOPO5.topo_sm_interp = topo_sm_interp;


%% pcolor_maker
fig_n = 1;
rowcols = [2 3];
rowcols_size = [8 5]; % cm
margs = [1 1 1 1]; % cm
gaps = [1 1]; % cm

cmaps_levels = 32;

for sp = 1 : rowcols(1)*rowcols(2)
    cmaps{sp} = flipud(othercolor('Paired12', cmaps_levels));
    minmax{sp} = [-2000 0];
    x{sp} = aus8_coor.lon;
    y{sp} = aus8_coor.lat;

end

x{1} = ETOPO5.lon;
y{1} = ETOPO5.lat;
x{4} = ETOPO5.lon;
y{4} = ETOPO5.lat;

topo_orig(topo_orig==0)=NaN;
topo_orig_interp(topo_orig_interp==0)=NaN;
topo_sm_interp(topo_sm_interp==0)=NaN;

data{1} = topo_orig;
data{2} = topo_orig_interp;
data{3} = topo_sm_interp;
data{4} = topo_orig;
data{5} = topo_orig_interp;
data{6} = topo_sm_interp;

titles{1} = ['ETOPO5 original bathym'];
titles{2} = ['ETOPO5 original bathym interp'];
titles{3} = ['ETOPO5 smoothed bathym interp $\Delta xy=' ...
    num2str(ETOPO5.sm_window) '^{\circ}$'];
titles{4} = ['GAB'];
titles{5} = ['GAB'];
titles{6} = ['GAB'];

axis_setup{1} = ...
    [aus8_coor.lon(1) aus8_coor.lon(end) ...
    aus8_coor.lat(end) aus8_coor.lat(1)];
axis_setup{2} = ...
    [aus8_coor.lon(1) aus8_coor.lon(end) ...
    aus8_coor.lat(end) aus8_coor.lat(1)];
axis_setup{3} = ...
    [aus8_coor.lon(1) aus8_coor.lon(end) ...
    aus8_coor.lat(end) aus8_coor.lat(1)];
axis_setup{4} = ...
    [125 135 -36 -32];
axis_setup{5} = ...
    [125 135 -36 -32];
axis_setup{6} = ...
    [125 135 -36 -32];

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
save([data_path 'SACS_data/ETOPO5'], 'ETOPO5')
disp('ETOPO5 DONE')

