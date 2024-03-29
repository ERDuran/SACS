%% get KDS75
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])

Months = aus8_coor.Months;


%%
for t = 1 : 12
    [data.(Months{t}), ~, kds_att] = ...
        nc2mat([data_path ...
        'KDS75/KDS75_ncra_y103to109_' Months{t} '.nc'], ...
        'ALL');
    
    disp([Months{t} ' OK!'])
end


% fix axis
KDS75_thet.lat = flipud(data.Jan.yt_ocean);
KDS75_thet.lon = (data.Jan.xt_ocean + 360)';
KDS75_thet.depth = -data.Jan.st_ocean;


%%
KDS75_thet.mean = 0;
KDS75_salt.mean = 0;

for t = 1 : 12
    KDS75_thet.(Months{t}) = ...
        flipud(permute(data.(Months{t}).temp, [2, 1, 3]));
    KDS75_salt.(Months{t}) = ...
        flipud(permute(data.(Months{t}).salt, [2, 1, 3]));
end
mean_thet = KDS75_thet.Jan;
mean_salt = KDS75_salt.Jan;
for t = 2 : 12
    mean_thet = mean_thet + KDS75_thet.(Months{t});
    mean_salt = mean_salt + KDS75_salt.(Months{t});
end
KDS75_thet.mean = mean_thet/12;
KDS75_salt.mean = mean_salt/12;


%% mean
fig_n = 1;
rowcols = [3 2];
rowcols_size = [8 5]; % cm
margs = [1 1 1 1]; % cm
gaps = [1 1]; % cm

for sp = 1 : rowcols(1)*rowcols(2)
    axis_setup{sp} = ...
        [KDS75_thet.lon(1) KDS75_thet.lon(end) ...
        KDS75_thet.lat(end) KDS75_thet.lat(1)];
    x{sp} = KDS75_thet.lon;
    y{sp} = KDS75_thet.lat;
end
    
z = {0 0 -400 -400 -1000 -1000};
for t = 1 : length(z)
    z_ind{t} = find_nearest(KDS75_thet.depth,z{t});
end

data = {...
    KDS75_thet.mean(:,:,z_ind{1}), ...
    KDS75_salt.mean(:,:,z_ind{2}), ...
    KDS75_thet.mean(:,:,z_ind{3}), ...
    KDS75_salt.mean(:,:,z_ind{4}), ...
    KDS75_thet.mean(:,:,z_ind{5}), ...
    KDS75_salt.mean(:,:,z_ind{6})};
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
    ['KDS75 mean $thet$ $z=' ...
    num2str(KDS75_thet.depth(z_ind{1})) '$ ($^{\circ}C$)'], ...
    ['KDS75 mean $salt$ $z=' ...
    num2str(KDS75_thet.depth(z_ind{2})) '$ ($psu$)'], ...
    ['KDS75 mean $thet$ $z=' ...
    num2str(KDS75_thet.depth(z_ind{3})) '$ ($^{\circ}C$)'], ...
    ['KDS75 mean $salt$ $z=' ...
    num2str(KDS75_thet.depth(z_ind{4})) '$ ($psu$)'], ...
    ['KDS75 mean $thet$ $z=' ...
    num2str(KDS75_thet.depth(z_ind{5})) '$ ($^{\circ}C$)'], ...
    ['KDS75 mean $salt$ $z=' ...
    num2str(KDS75_thet.depth(z_ind{6})) '$ ($psu$)']};
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
for t = 1 : length(z)
    z_ind{t} = find_nearest(KDS75_thet.depth,z{t});
end

month_ind = {'Jan', 'Apr', 'Jul', 'Oct', 'Jan', 'Apr', 'Jul', 'Oct', ...
    'Jan', 'Apr', 'Jul', 'Oct',};
cmaps_levels = 12;

for sp = 1 : rowcols(1)*rowcols(2)
    x{sp} = KDS75_thet.lon;
    y{sp} = KDS75_thet.lat;
    axis_setup{sp} = ...
        [KDS75_thet.lon(1) KDS75_thet.lon(end) ...
        KDS75_thet.lat(end) KDS75_thet.lat(1)];
    data{sp} = KDS75_thet.(month_ind{sp})(:,:,z_ind{sp});
    titles{sp} = ['KDS75 ' month_ind{sp} ...
        ' $thet$ $z=' ...
        num2str(KDS75_thet.depth(z_ind{sp})) '$ ($^{\circ}C$)'];
    cmaps{sp} = flipud(othercolor('RdYlBu8', cmaps_levels));
end

minmax = {...
    [8 20], ...
    [8 20], ...
    [8 20], ...
    [8 20], ...
    [min(min(data{5}))+2 max(max(data{5}))-2], ...
    [min(min(data{5}))+2 max(max(data{5}))-2], ...
    [min(min(data{5}))+2 max(max(data{5}))-2], ...
    [min(min(data{5}))+2 max(max(data{5}))-2], ...
    [min(min(data{9})) max(max(data{9}))-0.5], ...
    [min(min(data{9})) max(max(data{9}))-0.5], ...
    [min(min(data{9})) max(max(data{9}))-0.5], ...
    [min(min(data{9})) max(max(data{9}))-0.5]};

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
for t = 1 : length(z)
    z_ind{t} = find_nearest(KDS75_thet.depth,z{t});
end

month_ind = {'Jan', 'Apr', 'Jul', 'Oct', 'Jan', 'Apr', 'Jul', 'Oct', ...
    'Jan', 'Apr', 'Jul', 'Oct',};
cmaps_levels = 12;

for sp = 1 : rowcols(1)*rowcols(2)
    x{sp} = KDS75_thet.lon;
    y{sp} = KDS75_thet.lat;
    axis_setup{sp} = ...
        [KDS75_thet.lon(1) KDS75_thet.lon(end) ...
        KDS75_thet.lat(end) KDS75_thet.lat(1)];
    data{sp} = KDS75_salt.(month_ind{sp})(:,:,z_ind{sp});
    titles{sp} = ['KDS75 ' month_ind{sp} ...
        ' $salt$ $z=' ...
        num2str(KDS75_thet.depth(z_ind{sp})) '$ ($psu$)'];
    cmaps{sp} = flipud(othercolor('PuOr8', cmaps_levels));
end

minmax = {...
    [34 36], ...
    [34 36], ...
    [34 36], ...
    [34 36], ...
    [min(min(data{5}))+0.2 max(max(data{5}))-0.2], ...
    [min(min(data{5}))+0.2 max(max(data{5}))-0.2], ...
    [min(min(data{5}))+0.2 max(max(data{5}))-0.2], ...
    [min(min(data{5}))+0.2 max(max(data{5}))-0.2], ...
    [min(min(data{9})) max(max(data{9}))], ...
    [min(min(data{9})) max(max(data{9}))], ...
    [min(min(data{9})) max(max(data{9}))], ...
    [min(min(data{9})) max(max(data{9}))]};

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
save([data_path 'SACS_data/KDS75_thet'], 'KDS75_thet')
save([data_path 'SACS_data/KDS75_salt'], 'KDS75_salt')
disp('KDS75_thet KDS75_salt DONE')

