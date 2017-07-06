%% make Smith and Sandwell topography mask
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SmSan02/smsan_indian'])

lat = flipud(lat);
lon = lon';
topog = flipud(topog);


%%
lat_range = [-30 -50];
lon_range = [108 153];

lat_aus8 = (lat_range(1) : -1/8 : lat_range(2))';
lon_aus8 = lon_range(1) : 1/8 : lon_range(2);

% get latitudes within the range
lat_ind = find(lat <= lat_range(1)+1 & lat >= lat_range(2)-1);
% flip latitude indexes upside down
SmSan02.lat = lat(lat_ind);

% get longitudes within the range
lon_ind = find(lon >= lon_range(1)-1 & lon <= lon_range(2)+1);
% transpose lon to get it as a row vector
SmSan02.lon = lon(lon_ind)';

topo_orig = topog(lat_ind, lon_ind);

% non-smoothed
topo_orig(topo_orig>=0) = 0;
% non-smoothed
topo_orig_interp = interp2(...
    SmSan02.lon, SmSan02.lat, topo_orig, ...
    lon_aus8, lat_aus8, 'linear');

% smoothed
topo_orig(topo_orig==0) = NaN;
% we need to smooth out ETOPO-5 because there are some weird discrepancies
SmSan02.sm_factor = 2;
SmSan02.sm_window = (2*SmSan02.sm_factor+1)*1/30;
topo_sm = smooth2a(topo_orig,SmSan02.sm_factor);
% smoothed
topo_sm(isnan(topo_sm)) = 0;
topo_sm_interp = interp2(...
    SmSan02.lon, SmSan02.lat, topo_sm, lon_aus8, lat_aus8, 'linear');

topo_orig(topo_orig>=0) = NaN;
topo_orig_interp(topo_orig_interp>=0) = NaN;
topo_sm_interp(topo_sm_interp>=0) = NaN;


%% pcolor_maker
fig_n = 1;
rowcols = [2 3];
rowcols_size = [16 10]; % cm
margs = [1 1 1 1]; % cm
gaps = [1 1]; % cm

cmaps_levels = 10;

minmax{1} = [-2000 0];
minmax{2} = [-2000 0];
minmax{3} = [-2000 0];
minmax{4} = [-100 0];
minmax{5} = [-100 0];
minmax{6} = [-100 0];

for sp = 1 : rowcols(1)*rowcols(2)
    cmaps{sp} = flipud(othercolor('Paired12', cmaps_levels));
    x{sp} = aus8_coor.lon;
    y{sp} = aus8_coor.lat;

end

x{1} = SmSan02.lon;
y{1} = SmSan02.lat;
x{4} = SmSan02.lon;
y{4} = SmSan02.lat;

data{1} = topo_orig;
data{2} = topo_orig_interp;
data{3} = topo_sm_interp;
data{4} = topo_orig;
data{5} = topo_orig_interp;
data{6} = topo_sm_interp;

titles{1} = ['Smith and Sandwell 2002 original bathym'];
titles{2} = ['Smith and Sandwell 2002 original bathym interp'];
titles{3} = ['Smith and Sandwell 2002 ' ...
    'smoothed bathym interp $\Delta xy=' ...
    num2str(SmSan02.sm_window) '^{\circ}$'];
titles{4} = ['Bass Strait'];
titles{5} = ['Bass Strait'];
titles{6} = ['Bass Strait'];

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
    [140 150 -42 -38];
axis_setup{5} = ...
    [140 150 -42 -38];
axis_setup{6} = ...
    [140 150 -42 -38];

font_size = 9;
fig_color = [0 0 0];

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
%
topo_orig(isnan(topo_orig)) = 0;
topo_orig_interp(isnan(topo_orig_interp)) = 0;
topo_sm(isnan(topo_sm)) = 0;
topo_sm_interp(isnan(topo_sm_interp)) = 0;

%
SmSan02.topo_orig = topo_orig;
SmSan02.topo_orig_interp = topo_orig_interp;
SmSan02.topo_sm = topo_sm;
SmSan02.topo_sm_interp = topo_sm_interp;

save([data_path 'SACS_data/SmSan02'], 'SmSan02')
disp('SmSan02 DONE')

