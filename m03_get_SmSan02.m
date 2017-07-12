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

% 1) original topography
topo_orig = topog(lat_ind, lon_ind);
land_idx = topo_orig>=0;

% 2) directly interpolated topography
topo_orig(land_idx) = 0;
topo_orig_interp = interp2(...
    SmSan02.lon, SmSan02.lat, topo_orig, ...
    lon_aus8, lat_aus8, 'linear');
topo_orig_interp(topo_orig_interp==0) = NaN;

% 3) bin averaged interpolated topography
thrh_binavg = 8;
topo_orig(land_idx) = NaN;
topo_binavg = NaN(length(lat_aus8), length(lon_aus8));
for m = 1 : length(lat_aus8)
    lat_ind = find(SmSan02.lat <= lat_aus8(m)+1/16 & ...
        SmSan02.lat >= lat_aus8(m)-1/16);
    for n = 1 : length(lon_aus8)
        lon_ind = find(SmSan02.lon >= lon_aus8(n)-1/16 & ...
            SmSan02.lon <= lon_aus8(n)+1/16);
        binned_topo = topo_orig(lat_ind,lon_ind);
        binned_topo_finite = length(find(isfinite(binned_topo)));
        if binned_topo_finite > thrh_binavg 
            topo_binavg(m,n) = nanmean(nanmean(binned_topo));
        end
    end
    fprintf('lat = %5.3f\n', lat_aus8(m))
end


%% pcolor_maker
fig_n = 1;
rowcols = [3 3];
rowcols_size = [12 7.5]; % cm
margs = [1.5 1.5 1.5 1.5]; % cm
gaps = [1.5 1.5]; % cm

cmaps_levels = 10;

minmax{1} = [-2000 0];
minmax{2} = [-2000 0];
minmax{3} = [-2000 0];
minmax{4} = [-200 0];
minmax{5} = [-200 0];
minmax{6} = [-200 0];
minmax{7} = [-100 0];
minmax{8} = [-100 0];
minmax{9} = [-100 0];

for sp = 1 : rowcols(1)*rowcols(2)
    cmaps{sp} = flipud(othercolor('Paired12', cmaps_levels));
    x{sp} = aus8_coor.lon;
    y{sp} = aus8_coor.lat;
end

x{1} = SmSan02.lon;
y{1} = SmSan02.lat;
x{4} = SmSan02.lon;
y{4} = SmSan02.lat;
x{7} = SmSan02.lon;
y{7} = SmSan02.lat;

data{1} = topo_orig;
data{2} = topo_orig_interp;
data{3} = topo_binavg;
data{4} = topo_orig;
data{5} = topo_orig_interp;
data{6} = topo_binavg;
data{7} = topo_orig;
data{8} = topo_orig_interp;
data{9} = topo_binavg;

titles{1} = ['Smith and Sandwell 2002 original bathym'];
titles{2} = ['Smith and Sandwell 2002 original bathym interp'];
titles{3} = ['Smith and Sandwell 2002 bin averaged into aus8. ' ...
    'Threshold $>' num2str(thrh_binavg) '/16$'];
titles{4} = ['Great Australian Bight'];
titles{5} = ['Great Australian Bight'];
titles{6} = ['Great Australian Bight'];
titles{7} = ['Bass Strait'];
titles{8} = ['Bass Strait'];
titles{9} = ['Bass Strait'];

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
    [123 140 -37 -31];
axis_setup{5} = ...
    [123 140 -37 -31];
axis_setup{6} = ...
    [123 140 -37 -31];
axis_setup{7} = ...
    [140 150 -42 -38];
axis_setup{8} = ...
    [140 150 -42 -38];
axis_setup{9} = ...
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
topo_binavg(isnan(topo_binavg)) = 0;

%
SmSan02.topo_orig = topo_orig;
SmSan02.topo_orig_interp = topo_orig_interp;
SmSan02.topo_binavg = topo_binavg;

save([data_path 'SACS_data/SmSan02'], 'SmSan02')
disp('SmSan02 DONE')

