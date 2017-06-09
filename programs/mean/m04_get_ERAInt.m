%% get wind data
clearvars('-except', '*_path')

% wind
[wind_data, wind_gen_att, wind_att] = ...
    nc2mat([data_path ...
    'ERA_INTERIM/mth_av_turb_surf_stress_20032012/' ...
    '_grib2netcdf-atls12-a562cefde8a29a7288fa0b8b7f9413f7-qJ0lae.nc'], ...
    'ALL');

load([data_path 'SACS_data/ETOPO5'])


%%
ERAInt.lat = double(wind_data.latitude);
ERAInt.lon = double(wind_data.longitude)';
ERAInt.time_num = double(wind_data.time);
ERAInt.time_vec = datestr(datenum(1900,1,1) + ERAInt.time_num/24);

%
inss = permute(wind_data.inss, [2,1,3]);
ERAInt.Tau_y.mean = mean(inss,3);

%
iews = permute(wind_data.iews, [2,1,3]);
ERAInt.Tau_x.mean = mean(iews,3);

%
month_names = {...
    'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
month_vector = 1 : 12 : length(ERAInt.time_num);

for ll = 1 : 12
    ERAInt.Tau_y.(month_names{ll}) = ...
        mean(inss(:,:,month_vector+ll-1),3);
    ERAInt.Tau_x.(month_names{ll}) = ...
        mean(iews(:,:,month_vector+ll-1),3);
end


%% figure set-up
close all
font_size = 10;
fig_n = 1;
fig = figure(fig_n);
rowN = 3; colN = 1;
set(gcf,'units','centimeters','position',[0 0 colN*20 rowN*7]);

col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = 0.02; % gap width between subplots
gap_h = 0.04; % gap height between subplots
marg_b = 0.04; % bottom margin
marg_t = 0.04; % top margin
marg_l = 0.04; % left margin
marg_r = 0.02; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length
h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

% 1)
sp = 1;
axes(h_axes_sp(sp))
data1 = ERAInt.Tau_x.mean(:,:);
data1(ETOPO5.topo_sm_interp==0) = NaN;
data_min = -0.3;
data_max = 0.3;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdBu10', levels)));
pcolor(...
    ERAInt.lon, ...
    ERAInt.lat, ...
    data1)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['mean $\tau_{x}$ $(m^{2}/s)$'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 2)
sp = 2;
axes(h_axes_sp(sp))
data2 = ERAInt.Tau_y.mean(:,:);
data2(ETOPO5.topo_sm_interp==0) = NaN;
data_min = -0.1;
data_max = 0.1;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdBu10', levels)));
pcolor(...
    ERAInt.lon, ...
    ERAInt.lat, ...
    data2)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['mean $\tau_{y}$ $(m^{2}/s)$'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 3)
sp = 3;
axes(h_axes_sp(sp))
data = sqrt(data1.^2 + data2.^2);
data(ETOPO5.topo_sm_interp==0) = NaN;
data_min = 0;
data_max = 0.25;
levels = 12;
colormap(h_axes_sp(sp), othercolor('BuPu9', levels));
pcolor(...
    ERAInt.lon, ...
    ERAInt.lat, ...
    data)
shading interp
hold on
qn = 9; qs = 1;
quiver(...
    ERAInt.lon(1 : qn : end), ...
    ERAInt.lat(1 : qn : end), ...
    data1(1 : qn : end, 1 : qn : end), ...
    data2(1 : qn : end, 1 : qn : end), ...
    qs, 'k')
caxis([data_min data_max]);
cbar = colorbar;
title(['mean $\tau$ $(m^{2}/s)$'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% Save
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
close all
font_size = 10;
fig_n = 2;
fig = figure(fig_n);
rowN = 3; colN = 1;
set(gcf,'units','centimeters','position',[0 0 colN*20 rowN*7]);

col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = 0.02; % gap width between subplots
gap_h = 0.04; % gap height between subplots
marg_b = 0.04; % bottom margin
marg_t = 0.04; % top margin
marg_l = 0.04; % left margin
marg_r = 0.02; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length
h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

% 1)
sp = 1;
axes(h_axes_sp(sp))
data1 = ERAInt.Tau_x.Jan(:,:);
data1(ETOPO5.topo_sm_interp==0) = NaN;
data_min = -0.3;
data_max = 0.3;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdBu10', levels)));
pcolor(...
    ERAInt.lon, ...
    ERAInt.lat, ...
    data1)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['Jan $\tau_{x}$ $(m^{2}/s)$'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 2)
sp = 2;
axes(h_axes_sp(sp))
data2 = ERAInt.Tau_y.Jan(:,:);
data2(ETOPO5.topo_sm_interp==0) = NaN;
data_min = -0.1;
data_max = 0.1;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdBu10', levels)));
pcolor(...
    ERAInt.lon, ...
    ERAInt.lat, ...
    data2)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['Jan $\tau_{x}$ $(m^{2}/s)$'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 3)
sp = 3;
axes(h_axes_sp(sp))
data = sqrt(data1.^2 + data2.^2);
data(ETOPO5.topo_sm_interp==0) = NaN;
data_min = 0;
data_max = 0.25;
levels = 12;
colormap(h_axes_sp(sp), othercolor('BuPu9', levels));
pcolor(...
    ERAInt.lon, ...
    ERAInt.lat, ...
    data)
shading interp
hold on
qn = 9; qs = 1;
quiver(...
    ERAInt.lon(1 : qn : end), ...
    ERAInt.lat(1 : qn : end), ...
    data1(1 : qn : end, 1 : qn : end), ...
    data2(1 : qn : end, 1 : qn : end), ...
    qs, 'k')
caxis([data_min data_max]);
cbar = colorbar;
title(['Jan $\tau_{x}$ $(m^{2}/s)$'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% Save
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
close all
font_size = 10;
fig_n = 3;
fig = figure(fig_n);
rowN = 3; colN = 1;
set(gcf,'units','centimeters','position',[0 0 colN*20 rowN*7]);

col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = 0.02; % gap width between subplots
gap_h = 0.04; % gap height between subplots
marg_b = 0.04; % bottom margin
marg_t = 0.04; % top margin
marg_l = 0.04; % left margin
marg_r = 0.02; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length
h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

% 1)
sp = 1;
axes(h_axes_sp(sp))
data1 = ERAInt.Tau_x.Jul(:,:);
data1(ETOPO5.topo_sm_interp==0) = NaN;
data_min = -0.3;
data_max = 0.3;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdBu10', levels)));
pcolor(...
    ERAInt.lon, ...
    ERAInt.lat, ...
    data1)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['Jul $\tau_{x}$ $(m^{2}/s)$'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 2)
sp = 2;
axes(h_axes_sp(sp))
data2 = ERAInt.Tau_y.Jul(:,:);
data2(ETOPO5.topo_sm_interp==0) = NaN;
data_min = -0.1;
data_max = 0.1;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdBu10', levels)));
pcolor(...
    ERAInt.lon, ...
    ERAInt.lat, ...
    data2)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['Jul $\tau_{x}$ $(m^{2}/s)$'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 3)
sp = 3;
axes(h_axes_sp(sp))
data = sqrt(data1.^2 + data2.^2);
data(ETOPO5.topo_sm_interp==0) = NaN;
data_min = 0;
data_max = 0.25;
levels = 12;
colormap(h_axes_sp(sp), othercolor('BuPu9', levels));
pcolor(...
    ERAInt.lon, ...
    ERAInt.lat, ...
    data)
shading interp
hold on
qn = 9; qs = 1;
quiver(...
    ERAInt.lon(1 : qn : end), ...
    ERAInt.lat(1 : qn : end), ...
    data1(1 : qn : end, 1 : qn : end), ...
    data2(1 : qn : end, 1 : qn : end), ...
    qs, 'k')
caxis([data_min data_max]);
cbar = colorbar;
title(['Jul $\tau_{x}$ $(m^{2}/s)$'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% Save
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

