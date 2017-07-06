%% get KDS75
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])

Months = aus8_coor.Months;
lat = aus8_coor.lat;
lon = aus8_coor.lon;


%%
for t = 1 : 12
    [KDdata.(Months{t}), ~, kds_att] = ...
        nc2mat([data_path ...
        'KDS75/KDS75_ncra_y103to109_' Months{t} '.nc'], ...
        'ALL');
    
    disp([Months{t} ' OK!'])
end

yu = flipud(KDdata.Jan.yu_ocean);
xu = (KDdata.Jan.xu_ocean + 360)';


%%
KDau_tau.tau_x.mean = 0;
    KDau_tau.tau_y.mean = 0;

for t = 1 : 12
    KDau_tau.tau_x.(Months{t}) = interp2(xu, yu, ...
        flipud(permute(KDdata.(Months{t}).tau_x, [2, 1, 3])), lon, lat);
    
    KDau_tau.tau_y.(Months{t}) = interp2(xu, yu, ...
        flipud(permute(KDdata.(Months{t}).tau_y, [2, 1, 3])), lon, lat);
end

mean_tau_x = KDau_tau.tau_x.Jan;
mean_tau_y = KDau_tau.tau_y.Jan;
for t = 2 : 12
    mean_tau_x = mean_tau_x + KDau_tau.tau_x.(Months{t});
    mean_tau_y = mean_tau_y + KDau_tau.tau_y.(Months{t});
end
KDau_tau.tau_x.mean = mean_tau_x/12;
KDau_tau.tau_y.mean = mean_tau_y/12;


%% figure set-up
fig_n = 1;
rowcols = [3 1];
rowcols_size = [8 5]; % cm
margs = [1 1 1 1]; % cm
gaps = [1 1]; % cm

z_ind = {0 0 0 0 -400 -400 -400 -400 -1000 -1000 -1000 -1000};
month_ind = {'Jan', 'Apr', 'Jul', 'Oct', 'Jan', 'Apr', 'Jul', 'Oct', ...
    'Jan', 'Apr', 'Jul', 'Oct',};
cmaps_levels = 12;

for sp = 1 : rowcols(1)*rowcols(2)
    axis_setup{sp} = ...
        [aus8_coor.lon(1) aus8_coor.lon(end) ...
        aus8_coor.lat(end) aus8_coor.lat(1)];
    x{sp} = aus8_coor.lon;
    y{sp} = aus8_coor.lat;
    cmaps{sp} = flipud(othercolor('RdBu8', cmaps_levels));
    u{sp} = NaN(size(KDau_tau.tau_x.mean));
    v{sp} = NaN(size(KDau_tau.tau_x.mean));
end

u{3} = KDau_tau.tau_x.mean;
v{3} = KDau_tau.tau_y.mean;
n = 12;
s = 2;

data{1} = KDau_tau.tau_x.mean;
data{2} = KDau_tau.tau_y.mean;
data{3} = sqrt(data{1}.^2 + data{2}.^2);

titles{1} = ['KDS75 mean $\tau_{x}$ $(m^{2}/s)$'];
titles{2} = ['KDS75 mean $\tau_{y}$ $(m^{2}/s)$'];
titles{3} = ['KDS75 mean $\tau$ $(m^{2}/s)$'];

minmax = {...
    [-0.3 0.3], ...
    [-0.1 0.1], ...
    [0 0.25]};

cmaps{3} = othercolor('BuPu9', cmaps_levels);

font_size = 9;
fig_color = [0.7 0.7 0.7];

fig = quiver_maker(...
    fig_n, rowcols, rowcols_size, margs, gaps, ...
    x, y, data, u, v, n, s, axis_setup, minmax, cmaps, ...
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


%% figure set-up
fig_n = 2;
rowcols = [3 4];
rowcols_size = [8 5]; % cm
margs = [1 1 1 1]; % cm
gaps = [1 1]; % cm

z_ind = {0 0 0 0 -400 -400 -400 -400 -1000 -1000 -1000 -1000};
month_ind = {'Jan', 'Apr', 'Jul', 'Oct', 'Jan', 'Apr', 'Jul', 'Oct', ...
    'Jan', 'Apr', 'Jul', 'Oct',};
cmaps_levels = 12;

for sp = 1 : rowcols(1)*rowcols(2)
    axis_setup{sp} = ...
        [aus8_coor.lon(1) aus8_coor.lon(end) ...
        aus8_coor.lat(end) aus8_coor.lat(1)];
    x{sp} = aus8_coor.lon;
    y{sp} = aus8_coor.lat;
    cmaps{sp} = flipud(othercolor('RdBu8', cmaps_levels));
    u{sp} = NaN(size(KDau_tau.tau_x.mean));
    v{sp} = NaN(size(KDau_tau.tau_x.mean));
end

u{9} = KDau_tau.tau_x.Jan;
v{9} = KDau_tau.tau_y.Jan;
u{10} = KDau_tau.tau_x.Apr;
v{10} = KDau_tau.tau_y.Apr;
u{11} = KDau_tau.tau_x.Jul;
v{11} = KDau_tau.tau_y.Jul;
u{12} = KDau_tau.tau_x.Oct;
v{12} = KDau_tau.tau_y.Oct;
n = 12;
s = 2;


data{1} = KDau_tau.tau_x.Jan;
data{5} = KDau_tau.tau_y.Jan;
data{9} = sqrt(data{1}.^2 + data{5}.^2);

data{2} = KDau_tau.tau_x.Apr;
data{6} = KDau_tau.tau_y.Apr;
data{10} = sqrt(data{2}.^2 + data{6}.^2);

data{3} = KDau_tau.tau_x.Jul;
data{7} = KDau_tau.tau_y.Jul;
data{11} = sqrt(data{3}.^2 + data{7}.^2);

data{4} = KDau_tau.tau_x.Oct;
data{8} = KDau_tau.tau_y.Oct;
data{12} = sqrt(data{4}.^2 + data{8}.^2);

titles{1} = ['KDS75 Jan $\tau_{x}$ $(m^{2}/s)$ (03-12)'];
titles{5} = ['KDS75 Jan $\tau_{y}$ $(m^{2}/s)$ (03-12)'];
titles{9} = ['KDS75 Jan $\tau$ $(m^{2}/s)$ (03-12)'];
titles{2} = ['KDS75 Apr $\tau_{x}$ $(m^{2}/s)$ (03-12)'];
titles{6} = ['KDS75 Apr $\tau_{y}$ $(m^{2}/s)$ (03-12)'];
titles{10} = ['KDS75 Apr $\tau$ $(m^{2}/s)$ (03-12)'];
titles{3} = ['KDS75 Jul $\tau_{x}$ $(m^{2}/s)$ (03-12)'];
titles{7} = ['KDS75 Jul $\tau_{y}$ $(m^{2}/s)$ (03-12)'];
titles{11} = ['KDS75 Jul $\tau$ $(m^{2}/s)$ (03-12)'];
titles{4} = ['KDS75 Oct $\tau_{x}$ $(m^{2}/s)$ (03-12)'];
titles{8} = ['KDS75 Oct $\tau_{y}$ $(m^{2}/s)$ (03-12)'];
titles{12} = ['KDS75 Oct $\tau$ $(m^{2}/s)$ (03-12)'];

minmax = {...
    [-0.3 0.3], ...
    [-0.3 0.3], ...
    [-0.3 0.3], ...
    [-0.3 0.3], ...
    [-0.1 0.1], ...
    [-0.1 0.1], ...
    [-0.1 0.1], ...
    [-0.1 0.1], ...
    [0 0.25], ...
    [0 0.25], ...
    [0 0.25], ...
    [0 0.25]};

cmaps{9} = othercolor('BuPu9', cmaps_levels);
cmaps{10} = othercolor('BuPu9', cmaps_levels);
cmaps{11} = othercolor('BuPu9', cmaps_levels);
cmaps{12} = othercolor('BuPu9', cmaps_levels);

font_size = 9;
fig_color = [0.7 0.7 0.7];

fig = quiver_maker(...
    fig_n, rowcols, rowcols_size, margs, gaps, ...
    x, y, data, u, v, n, s, axis_setup, minmax, cmaps, ...
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
save([data_path 'SACS_data/KDau_tau'], 'KDau_tau')
disp('KDau_tau DONE')


