%% Get aus8 data at selected lon lat depth
clearvars('-except', '*_path')

% temperature
[temp_data, temp_gen_att, temp_att] = ...
    nc2mat([data_path ...
    'CARS/aus8/temperature_aus8_2013.nc'], ...
    'ALL');

% salinity
[salt_data, salt_gen_att, salt_att] = ...
    nc2mat([data_path ...
    'CARS/aus8/salinity_aus8_2013.nc'], ...
    'ALL');

lat = double(temp_data.lat);
lon = double(temp_data.lon);
depth = -double(temp_data.depth);
depth_ann = -double(temp_data.depth_ann);
depth_semiann = -double(temp_data.depth_semiann);

temp.mean = temp_data.mean;
temp.nq = temp_data.nq;
temp.radius_q = temp_data.radius_q;
temp.an_cos = temp_data.an_cos;
temp.an_sin = temp_data.an_sin;
temp.sa_cos = temp_data.sa_cos;
temp.sa_sin = temp_data.sa_sin;

salt.mean = salt_data.mean;
salt.nq = salt_data.nq;
salt.radius_q = salt_data.radius_q;
salt.an_cos = salt_data.an_cos;
salt.an_sin = salt_data.an_sin;
salt.sa_cos = salt_data.sa_cos;
salt.sa_sin = salt_data.sa_sin;


%% Subset variables of aus8 from lon and lat window
% Define latitude and longitude ranges to zoom in
lat_range = [-30 -50];
lon_range = [108 153];
depth_max = -2000;

% get latitudes within the range
lat_ind = find(lat <= lat_range(1) & lat >= lat_range(2));
% flip latitude indexes upside down
lat_ind = flipud(lat_ind);

% get longitudes within the range
lon_ind = find(lon >= lon_range(1) & lon <= lon_range(2));

depth_ind = find(depth >= depth_max);

aus8_coor.lat = lat(lat_ind);
aus8_coor.lon = lon(lon_ind)';
aus8_coor.depth = depth(depth_ind);
save([data_path 'SACS_data/aus8_coor'], 'aus8_coor')

aus8.depth_ann = depth_ann;
aus8.depth_semiann = depth_semiann;

aus8.temp.mean = permute(temp.mean(lon_ind,lat_ind,depth_ind),[2,1,3]);
aus8.temp.nq = permute(temp.nq(lon_ind,lat_ind,depth_ind),[2,1,3]);
aus8.temp.radius_q = ...
    permute(temp.radius_q(lon_ind,lat_ind,depth_ind),[2,1,3]);
aus8.temp.an_cos = permute(temp.an_cos(lon_ind,lat_ind,:),[2,1,3]);
aus8.temp.an_sin = permute(temp.an_sin(lon_ind,lat_ind,:),[2,1,3]);
aus8.temp.sa_cos = permute(temp.sa_cos(lon_ind,lat_ind,:),[2,1,3]);
aus8.temp.sa_sin = permute(temp.sa_sin(lon_ind,lat_ind,:),[2,1,3]);

aus8.salt.mean = permute(salt.mean(lon_ind,lat_ind,depth_ind),[2,1,3]);
aus8.salt.nq = permute(salt.nq(lon_ind,lat_ind,depth_ind),[2,1,3]);
aus8.salt.radius_q = ...
    permute(salt.radius_q(lon_ind,lat_ind,depth_ind),[2,1,3]);
aus8.salt.an_cos = permute(salt.an_cos(lon_ind,lat_ind,:),[2,1,3]);
aus8.salt.an_sin = permute(salt.an_sin(lon_ind,lat_ind,:),[2,1,3]);
aus8.salt.sa_cos = permute(salt.sa_cos(lon_ind,lat_ind,:),[2,1,3]);
aus8.salt.sa_sin = permute(salt.sa_sin(lon_ind,lat_ind,:),[2,1,3]);


% make sure temp and salt have the same dry points
temp_nan = isnan(aus8.temp.mean);
salt_nan = isnan(aus8.salt.mean);
aus8.temp.mean(salt_nan) = NaN;
aus8.salt.mean(temp_nan) = NaN;


%% figure set-up
close all
font_size = 10;
fig_n = 1;
fig = figure(fig_n);
rowN = 3; colN = 2;
set(gcf,'units','centimeters','position',[0 0 colN*10 rowN*7]);

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
depth_lvl = 0;
depth_lvl_ind = depth == depth_lvl;
data = aus8.temp.mean(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu11', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.temp.mean z=' num2str(depth_lvl) ' (degc)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 2)
sp = 2;
axes(h_axes_sp(sp))
depth_lvl = 0;
depth_lvl_ind = depth == depth_lvl;
data = aus8.salt.mean(:,:,depth_lvl_ind);
data_min = 33.5;
data_max = 37;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('YlGnBu9', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.salt.mean z=' num2str(depth_lvl) ' (psu)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 3)
sp = 3;
axes(h_axes_sp(sp))
depth_lvl = -400;
depth_lvl_ind = depth == depth_lvl;
data = aus8.temp.mean(:,:,depth_lvl_ind);
data_min = 6;
data_max = 12;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu11', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.temp.mean z=' num2str(depth_lvl) ' (degc)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 4)
sp = 4;
axes(h_axes_sp(sp))
depth_lvl = -400;
depth_lvl_ind = depth == depth_lvl;
data = aus8.salt.mean(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = 35.2;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('YlGnBu9', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.salt.mean z=' num2str(depth_lvl) ' (psu)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 5)
sp = 5;
axes(h_axes_sp(sp))
depth_lvl = -1000;
depth_lvl_ind = depth == depth_lvl;
data = aus8.temp.mean(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu11', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.temp.mean z=' num2str(depth_lvl) ' (degC)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 6)
sp = 6;
axes(h_axes_sp(sp))
depth_lvl = -1000;
depth_lvl_ind = depth == depth_lvl;
data = aus8.salt.mean(:,:,depth_lvl_ind);
data_min = 34.25;
data_max = 34.55;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('YlGnBu9', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.salt.mean z=' num2str(depth_lvl) ' (psu)'])
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


%%
save([data_path 'SACS_data/aus8'], 'aus8')
disp('aus8 DONE')

