%% make anomaly
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
aus8_coor_names = fieldnames(aus8_coor);
a = aus8_coor;
load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/KDau_temp'])
load([data_path 'SACS_data/KDau_salt'])
load([data_path 'SACS_data/aus8_thet'])
load([data_path 'SACS_data/aus8_salt'])


%%
a.Months{13}='mean';

for t = 1 : length(a.Months)
    b.thet.(a.Months{t}) = ...
        aus8_thet.(a.Months{t}) - KDau_temp.(a.Months{t});
    b.salt.(a.Months{t}) = ...
        aus8_salt.(a.Months{t}) - KDau_salt.(a.Months{t});
    
end


%%
fig_n = 1;
rowcols = [3 2];
rowcols_size = [8 5]; % cm
margs = [1 1 1 1]; % cm
gaps = [1 1]; % cm

for sp = 1 : rowcols(1)*rowcols(2)
    axis_setup{sp} = ...
        [a.lon(1) a.lon(end) ...
        a.lat(end) a.lat(1)];
    x{sp} = a.lon;
    y{sp} = a.lat;
end
    
z = {0 0 -400 -400 -1000 -1000};
for p = 1 : length(z)
    z_ind{p} = find_nearest(a.depth,z{p});
end

data = {...
    b.thet.mean(:,:,z_ind{1}), ...
    b.salt.mean(:,:,z_ind{2}), ...
    b.thet.mean(:,:,z_ind{3}), ...
    b.salt.mean(:,:,z_ind{4}), ...
    b.thet.mean(:,:,z_ind{5}), ...
    b.salt.mean(:,:,z_ind{6})};
minmax = {...
    [-1.5 1.5], ...
    [-0.6 0.6], ...
    [-1.5 1.5], ...
    [-0.4 0.4], ...
    [-1.5 1.5], ...
    [-0.2 0.2]};
cmaps_levels = 12;
cmaps = {...
    flipud(othercolor('RdYlBu8', cmaps_levels)), ...
    flipud(othercolor('PuOr8', cmaps_levels)), ...
    flipud(othercolor('RdYlBu8', cmaps_levels)), ...
    flipud(othercolor('PuOr8', cmaps_levels)), ...
    flipud(othercolor('RdYlBu8', cmaps_levels)), ...
    flipud(othercolor('PuOr8', cmaps_levels))};
titles = {...
    ['aus8 - KDS75 mean $temp$ $z=' ...
    num2str(aus8_coor.depth(z_ind{1})) '$ ($^{\circ}C$)'], ...
    ['aus8 - KDS75 mean $salt$ $z=' ...
    num2str(aus8_coor.depth(z_ind{2})) '$ ($psu$)'], ...
    ['aus8 - KDS75 mean $temp$ $z=' ...
    num2str(aus8_coor.depth(z_ind{3})) '$ ($^{\circ}C$)'], ...
    ['aus8 - KDS75 mean $salt$ $z=' ...
    num2str(aus8_coor.depth(z_ind{4})) '$ ($psu$)'], ...
    ['aus8 - KDS75 mean $temp$ $z=' ...
    num2str(aus8_coor.depth(z_ind{5})) '$ ($^{\circ}C$)'], ...
    ['aus8 - KDS75 mean $salt$ $z=' ...
    num2str(aus8_coor.depth(z_ind{6})) '$ ($psu$)']};
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


%% seas var temp
fig_n = 2;
rowcols = [3 4];
rowcols_size = [8 5]; % cm
margs = [1 1 1 1]; % cm
gaps = [1 1]; % cm

z = {0 0 0 0 -400 -400 -400 -400 -1000 -1000 -1000 -1000};
for p = 1 : length(z)
    z_ind{p} = find_nearest(a.depth,z{p});
end

month_ind = {'Jan', 'Apr', 'Jul', 'Oct', 'Jan', 'Apr', 'Jul', 'Oct', ...
    'Jan', 'Apr', 'Jul', 'Oct',};
cmaps_levels = 12;

for sp = 1 : rowcols(1)*rowcols(2)
    x{sp} = a.lon;
    y{sp} = a.lat;
    axis_setup{sp} = ...
        [a.lon(1) a.lon(end) ...
        a.lat(end) a.lat(1)];
    data{sp} = b.thet.(month_ind{sp})(:,:,z_ind{sp});
    titles{sp} = ['aus8 - KDS75 ' month_ind{sp} ...
        ' $temp$ $z=' ...
        num2str(a.depth(z_ind{sp})) '$ ($^{\circ}C$)'];
    cmaps{sp} = flipud(othercolor('RdYlBu8', cmaps_levels));
end

minmax = {...
    [-1.5 1.5], ...
    [-1.5 1.5], ...
    [-1.5 1.5], ...
    [-1.5 1.5], ...
    [-1.5 1.5], ...
    [-1.5 1.5], ...
    [-1.5 1.5], ...
    [-1.5 1.5], ...
    [-1.5 1.5], ...
    [-1.5 1.5], ...
    [-1.5 1.5], ...
    [-1.5 1.5]};

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


%% seas var salt
fig_n = 3;
rowcols = [3 4];
rowcols_size = [8 5]; % cm
margs = [1 1 1 1]; % cm
gaps = [1 1]; % cm

z = {0 0 0 0 -400 -400 -400 -400 -1000 -1000 -1000 -1000};
for p = 1 : length(z)
    z_ind{p} = find_nearest(a.depth,z{p});
end

month_ind = {'Jan', 'Apr', 'Jul', 'Oct', 'Jan', 'Apr', 'Jul', 'Oct', ...
    'Jan', 'Apr', 'Jul', 'Oct',};
cmaps_levels = 12;

for sp = 1 : rowcols(1)*rowcols(2)
    x{sp} = a.lon;
    y{sp} = a.lat;
    axis_setup{sp} = ...
        [a.lon(1) a.lon(end) ...
        a.lat(end) a.lat(1)];
    data{sp} = b.salt.(month_ind{sp})(:,:,z_ind{sp});
    titles{sp} = ['aus8 - KDS75 ' month_ind{sp} ...
        ' $salt$ $z=' ...
        num2str(a.depth(z_ind{sp})) '$ ($psu$)'];
    cmaps{sp} = flipud(othercolor('PuOr8', cmaps_levels));
end

minmax = {...
    [-0.6 0.6], ...
    [-0.6 0.6], ...
    [-0.6 0.6], ...
    [-0.6 0.6], ...
    [-0.4 0.4], ...
    [-0.4 0.4], ...
    [-0.4 0.4], ...
    [-0.4 0.4], ...
    [-0.2 0.2], ...
    [-0.2 0.2], ...
    [-0.2 0.2], ...
    [-0.2 0.2]};

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
