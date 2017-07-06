%% get wind data
clearvars('-except', '*_path')

% wind
[wind_data, wind_gen_att, wind_att] = ...
    nc2mat([data_path ...
    'ERA_INTERIM/mth_av_turb_surf_stress_20032012/' ...
    '_grib2netcdf-atls12-a562cefde8a29a7288fa0b8b7f9413f7-qJ0lae.nc'], ...
    'ALL');

load([data_path 'SACS_data/SmSan02'])
load([data_path 'SACS_data/aus8_coor'])


%%
ERAInt.lat = double(wind_data.latitude);
ERAInt.lon = double(wind_data.longitude)';
ERAInt.time_num = double(wind_data.time);
ERAInt.time_vec = datestr(datenum(1900,1,1) + ERAInt.time_num/24);

%
inss = permute(wind_data.inss, [2,1,3]);
ERAInt.Tau_y.mean = mean(inss,3);
ERAInt.Tau_y.mean(SmSan02.topo_sm_interp==0) = NaN;
%
iews = permute(wind_data.iews, [2,1,3]);
ERAInt.Tau_x.mean = mean(iews,3);
ERAInt.Tau_x.mean(SmSan02.topo_sm_interp==0) = NaN;

%
Months = aus8_coor.Months;
month_vector = 1 : 12 : length(ERAInt.time_num);

for ll = 1 : 12
    ERAInt.Tau_y.(Months{ll}) = ...
        mean(inss(:,:,month_vector+ll-1),3);
    ERAInt.Tau_y.(Months{ll})(SmSan02.topo_sm_interp==0) = NaN;
        
    ERAInt.Tau_x.(Months{ll}) = ...
        mean(iews(:,:,month_vector+ll-1),3);
    ERAInt.Tau_x.(Months{ll})(SmSan02.topo_sm_interp==0) = NaN;
end


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
        [ERAInt.lon(1) ERAInt.lon(end) ...
        ERAInt.lat(end) ERAInt.lat(1)];
    x{sp} = ERAInt.lon;
    y{sp} = ERAInt.lat;
    cmaps{sp} = flipud(othercolor('RdBu8', cmaps_levels));
    u{sp} = NaN(size(ERAInt.Tau_x.mean));
    v{sp} = NaN(size(ERAInt.Tau_x.mean));
    nn{sp} = 12;
    s{sp} = 2;
end

u{3} = ERAInt.Tau_x.mean;
v{3} = ERAInt.Tau_y.mean;

data{1} = ERAInt.Tau_x.mean;
data{2} = ERAInt.Tau_y.mean;
data{3} = sqrt(data{1}.^2 + data{2}.^2);

titles{1} = ['ERA-Int mean $\tau_{x}$ $(m^{2}/s)$ (03-12)'];
titles{2} = ['ERA-Int mean $\tau_{y}$ $(m^{2}/s)$ (03-12)'];
titles{3} = ['ERA-Int mean $\tau$ $(m^{2}/s)$ (03-12)'];

minmax = {...
    [-0.3 0.3], ...
    [-0.1 0.1], ...
    [0 0.25]};

cmaps{3} = othercolor('BuPu9', cmaps_levels);

font_size = 9;
fig_color = [0.7 0.7 0.7];

fig = quiver_maker(...
    fig_n, rowcols, rowcols_size, margs, gaps, ...
    x, y, data, u, v, nn, s, axis_setup, minmax, cmaps, ...
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
        [ERAInt.lon(1) ERAInt.lon(end) ...
        ERAInt.lat(end) ERAInt.lat(1)];
    x{sp} = ERAInt.lon;
    y{sp} = ERAInt.lat;
    cmaps{sp} = flipud(othercolor('RdBu8', cmaps_levels));
    u{sp} = NaN(size(ERAInt.Tau_x.mean));
    v{sp} = NaN(size(ERAInt.Tau_x.mean));
    nn{sp} = 12;
    s{sp} = 2;
end

u{9} = ERAInt.Tau_x.Jan;
v{9} = ERAInt.Tau_y.Jan;
u{10} = ERAInt.Tau_x.Apr;
v{10} = ERAInt.Tau_y.Apr;
u{11} = ERAInt.Tau_x.Jul;
v{11} = ERAInt.Tau_y.Jul;
u{12} = ERAInt.Tau_x.Oct;
v{12} = ERAInt.Tau_y.Oct;


data{1} = ERAInt.Tau_x.Jan;
data{1}(SmSan02.topo_sm_interp==0) = NaN;
data{5} = ERAInt.Tau_y.Jan;
data{5}(SmSan02.topo_sm_interp==0) = NaN;
data{9} = sqrt(data{1}.^2 + data{5}.^2);

data{2} = ERAInt.Tau_x.Apr;
data{2}(SmSan02.topo_sm_interp==0) = NaN;
data{6} = ERAInt.Tau_y.Apr;
data{6}(SmSan02.topo_sm_interp==0) = NaN;
data{10} = sqrt(data{2}.^2 + data{6}.^2);

data{3} = ERAInt.Tau_x.Jul;
data{3}(SmSan02.topo_sm_interp==0) = NaN;
data{7} = ERAInt.Tau_y.Jul;
data{7}(SmSan02.topo_sm_interp==0) = NaN;
data{11} = sqrt(data{3}.^2 + data{7}.^2);

data{4} = ERAInt.Tau_x.Oct;
data{4}(SmSan02.topo_sm_interp==0) = NaN;
data{8} = ERAInt.Tau_y.Oct;
data{8}(SmSan02.topo_sm_interp==0) = NaN;
data{12} = sqrt(data{4}.^2 + data{8}.^2);

titles{1} = ['ERA-Int Jan $\tau_{x}$ $(m^{2}/s)$ (03-12)'];
titles{5} = ['ERA-Int Jan $\tau_{y}$ $(m^{2}/s)$ (03-12)'];
titles{9} = ['ERA-Int Jan $\tau$ $(m^{2}/s)$ (03-12)'];
titles{2} = ['ERA-Int Apr $\tau_{x}$ $(m^{2}/s)$ (03-12)'];
titles{6} = ['ERA-Int Apr $\tau_{y}$ $(m^{2}/s)$ (03-12)'];
titles{10} = ['ERA-Int Apr $\tau$ $(m^{2}/s)$ (03-12)'];
titles{3} = ['ERA-Int Jul $\tau_{x}$ $(m^{2}/s)$ (03-12)'];
titles{7} = ['ERA-Int Jul $\tau_{y}$ $(m^{2}/s)$ (03-12)'];
titles{11} = ['ERA-Int Jul $\tau$ $(m^{2}/s)$ (03-12)'];
titles{4} = ['ERA-Int Oct $\tau_{x}$ $(m^{2}/s)$ (03-12)'];
titles{8} = ['ERA-Int Oct $\tau_{y}$ $(m^{2}/s)$ (03-12)'];
titles{12} = ['ERA-Int Oct $\tau$ $(m^{2}/s)$ (03-12)'];

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
    x, y, data, u, v, nn, s, axis_setup, minmax, cmaps, ...
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


%% save
save([data_path 'SACS_data/ERAInt'], 'ERAInt')
disp('ERAInt DONE')

