%% get KDS75
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])

Months = aus8_coor.Months;
lat = aus8_coor.lat;
lon = aus8_coor.lon;
lat_u = lat(1:end-1) - 1/16;
lon_u = lon;
lat_v = lat;
lon_v = lon(1:end-1) + 1/16;


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
        flipud(permute(KDdata.(Months{t}).tau_x, [2, 1, 3])), ...
        lon_v, lat_v);
    
    KDau_tau.tau_y.(Months{t}) = interp2(xu, yu, ...
        flipud(permute(KDdata.(Months{t}).tau_y, [2, 1, 3])), ...
        lon_u, lat_u);
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
    cmaps{sp} = flipud(othercolor('RdBu8', cmaps_levels));
end

x{1} = lon_v;
y{1} = lat_v;
x{2} = lon_u;
y{2} = lat_u;
x{3} = lon;
y{3} = lat;

n = {12 12 12};
s = {2 2 2};

data{1} = KDau_tau.tau_x.mean;
data{2} = KDau_tau.tau_y.mean;   
tau_x_interp = interp2(lon_v, lat_v, data{1}, lon, lat);
tau_y_interp = interp2(lon_u, lat_u, data{2}, lon, lat); 
data{3} = sqrt(tau_x_interp.^2 + tau_y_interp.^2);

u{1} = NaN(size(data{1}));
v{1} = NaN(size(data{1}));
u{2} = NaN(size(data{2}));
v{2} = NaN(size(data{2}));
u{3} = tau_x_interp;
v{3} = tau_y_interp;

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

cmaps_chc = {'RdBu8', 'BuPu8'};
cmaps_ind = [1 1 1 1 1 1 1 1 2 2 2 2];
x_chc = {lon_v, lon_u, lon};
y_chc = {lat_v, lat_u, lat};
xy_ind = [1 1 1 1 2 2 2 2 3 3 3 3];
minmax_chc = {...
    [-0.3 0.3], ...
    [-0.1 0.1], ...
    [0 0.25]};

for sp = 1 : rowcols(1)*rowcols(2)
    axis_setup{sp} = ...
        [aus8_coor.lon(1) aus8_coor.lon(end) ...
        aus8_coor.lat(end) aus8_coor.lat(1)];
    x{sp} = x_chc{xy_ind(sp)};
    y{sp} = y_chc{xy_ind(sp)};
    cmaps{sp} = ...
        othercolor(cmaps_chc{cmaps_ind(sp)}, cmaps_levels);
    if cmaps_ind(sp) == 1
        cmaps{sp} = flipud(cmaps{sp});
    end
    minmax{sp} = minmax_chc{xy_ind(sp)};
    n{sp} = 12;
    s{sp} = 2;
end

for sp = 1 : rowcols(1)*rowcols(2)
    if cmaps_ind(sp) == 1
        if xy_ind(sp) == 1
            u{sp} = NaN(size(KDau_tau.tau_x.mean));
            v{sp} = NaN(size(KDau_tau.tau_x.mean));
            data{sp} = KDau_tau.tau_x.(month_ind{sp});
            titles{sp} = ...
                ['KDS75 ' month_ind{sp} ' $\tau_{x}$ $(m^{2}/s)$'];
        elseif xy_ind(sp) == 2
            u{sp} = NaN(size(KDau_tau.tau_y.mean));
            v{sp} = NaN(size(KDau_tau.tau_y.mean));
            data{sp} = KDau_tau.tau_y.(month_ind{sp});
            titles{sp} = ...
                ['KDS75 ' month_ind{sp} ' $\tau_{y}$ $(m^{2}/s)$'];
        end
    elseif cmaps_ind(sp) == 2
        u{sp} = interp2(lon_v,lat_v, ...
                KDau_tau.tau_x.(month_ind{sp}),lon,lat);
        v{sp} = interp2(lon_u,lat_u, ...
                KDau_tau.tau_y.(month_ind{sp}),lon,lat);
        data{sp} = sqrt(u{sp}.^2 + v{sp}.^2);
        titles{sp} = ...
            ['KDS75 ' month_ind{sp} ' $\tau$ $(m^{2}/s)$'];
    end
end

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


