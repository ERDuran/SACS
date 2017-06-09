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
data = ETOPO5.topo_orig(:,:);
data_min = -2000;
data_max = 0;
levels = 20;
colormap(h_axes_sp(sp), othercolor('Cat_12', levels));
pcolor(...
    ETOPO5.lon, ...
    ETOPO5.lat, ...
    data)
shading interp
axis([125 135 -36 -32])
caxis([data_min data_max]);
cbar = colorbar;
title(['ETOPO5.orig'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 2)
sp = 2;
axes(h_axes_sp(sp))
data = ETOPO5.topo_orig_interp(:,:);
data_min = -2000;
data_max = 0;
levels = 20;
colormap(h_axes_sp(sp), othercolor('Cat_12', levels));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
axis([125 135 -36 -32])
caxis([data_min data_max]);
cbar = colorbar;
title('ETOPO5.orig.interp')
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 3)
sp = 3;
axes(h_axes_sp(sp))
data = ETOPO5.topo_sm_interp(:,:);
data_min = -2000;
data_max = 0;
levels = 20;
colormap(h_axes_sp(sp), othercolor('Cat_12', levels));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
axis([125 135 -36 -32])
caxis([data_min data_max]);
cbar = colorbar;
title(['ETOPO5.sm.interp sm.window=' num2str(ETOPO5.sm_window) 'deg'])
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
save([data_path 'SACS_data/ETOPO5'], 'ETOPO5')
disp('ETOPO5 DONE')

