%% fig 1: map at the surface of T in colours S in contours
% and V in arrows
clc
path(pathdef)

% set up main directory
cd ~
% cd /home/z5100174/Desktop/MATLAB

% Add path to the data to be loaded
addpath(genpath('Dropbox/SACS_work'))
addpath(genpath(['/Users/earl/Dropbox/' ...
    'LeeuwinUndercurrent_HonoursProject/' ...
    'matlab/OFAM/ofam_out']))

clear

load aus8_TEOS10
load aus8_geostrophy

load ofam_mean


%% load stuff
lat_CARS = aus8_TEOS10.lat;
lon_CARS = aus8_TEOS10.lon;
pres = aus8_TEOS10.pres;
Theta_CARS = aus8_TEOS10.Theta.mean;
asal_CARS = aus8_TEOS10.asal.mean;

lat_u_CARS = aus8_geostrophy.u.lat_u;
lon_u_CARS = lon_CARS;
u_CARS = aus8_geostrophy.u.u_0_HH;

lat_v_CARS = lat_CARS;
lon_v_CARS = aus8_geostrophy.v.lon_v;
v_CARS = aus8_geostrophy.v.v_0_HH;


lat_OFAM = ofam_mean.lat;
lon_OFAM = ofam_mean.lon;
pres = ofam_mean.pressure;
Theta_OFAM = ofam_mean.cons_temperature;
asal_OFAM = ofam_mean.abs_salinity;

lat_u_OFAM = lat_OFAM;
lon_u_OFAM = lon_OFAM;
u_OFAM = ofam_mean.u_vel;

lat_v_OFAM = lat_OFAM;
lon_v_OFAM = lon_OFAM;
v_OFAM = ofam_mean.v_vel;


%% plot map
% 1) figure set-up
close all
font_size = 8;
fig1 = figure(1);
% magn = 3;
% set(gcf,'PaperType','A4', ...
%     'paperOrientation', 'portrait', ...
%     'paperunits','centimeters', ...
%     'PaperPosition',[0, 0, (21-4)/2, (29.7-4)]*magn);
set(gcf, 'units', 'normalized', 'outerposition', [0 0 0.7 0.95]);
rowN = 2; colN = 2;
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = - 0.04; % gap width between subplots
gap_h = 0.04; % gap height between subplots
marg_b = 0.03; % bottom margin
marg_t = 0.03; % top margin
marg_l = 0.02; % left margin
marg_r = 0.02; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length
h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);
colormap(flipud(othercolor('RdBu11')))

% 2) all data set-up
sp = 1;
axes(h_axes_sp(sp))
pres_now = 0;
pres_now_ind = pres == pres_now;
lon_hr = lon_CARS(1) : 0.1 : lon_CARS(end);
lat_hr = (lat_CARS(1) : -0.1 : lat_CARS(end))';
u_interp2_hr_now = ...
    interp2(lon_u_CARS, lat_u_CARS, ...
    u_CARS(:,:,pres_now_ind), lon_hr, lat_hr);
v_interp2_hr_now = ...
    interp2(lon_v_CARS, lat_v_CARS, ...
    v_CARS(:,:,pres_now_ind), lon_hr, lat_hr);
Theta_now = Theta_CARS(:,:,pres_now_ind);
asal_now = asal_CARS(:,:,pres_now_ind);

% 3) asal pcolor set-up
asal_contours = 34.7:0.1:36.5;
cmap_levels = length(asal_contours);
cmap = flipud(othercolor('RdBu11', cmap_levels-1));
colormap(h_axes_sp(sp), cmap);

% 4) plot asal pcolor
pcolor(...
    lon_CARS, ...
    lat_CARS, ...
    asal_now)
shading interp
hold on
caxis([asal_contours(1) asal_contours(end)]);
axis([110 152 -47 -31])
cbar = colorbar;
set(cbar, 'visible', 'off');

% 5) Theta contours set-up
Theta_contours = 10:1:21;

% 6) plot Theta labelled contours
[c, h] = contour(...
    lon_CARS, ...
    lat_CARS, ...
    Theta_now, ...
    Theta_contours, ...
    'k', ...
    'linewidth',0.5);
clabel(c,h, ...
    'fontsize', font_size-1, ...
    'fontweight', 'bold')

% 7) UV quiver set-up
[lon_mg, lat_mg]=meshgrid(lon_hr,lat_hr);
quiv_S = 5;
nn = 2;
u_interp2_hr_now(u_interp2_hr_now<0) = NaN;

% 8) plot directional U and V
quiver(...
    lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
    u_interp2_hr_now(1:nn:end, 1:nn:end), ...
    v_interp2_hr_now(1:nn:end, 1:nn:end), ...
    quiv_S, 'k');
% reference arrow goes here

% 9) title, grid, background and fonts
title('CARS surface')
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end


% 10) 2nd subplot: OFAM @ surface
sp = 2;
axes(h_axes_sp(sp))
pres_now = 0;
pres_now_ind = pres == pres_now;
u_now = u_OFAM(:,:,pres_now_ind);
v_now = v_OFAM(:,:,pres_now_ind);
Theta_now = Theta_OFAM(:,:,pres_now_ind);
asal_now = asal_OFAM(:,:,pres_now_ind);
asal_contours = 34.7:0.1:36.5;
cmap_levels = length(asal_contours);
cmap = flipud(othercolor('RdBu11', cmap_levels-1));
colormap(h_axes_sp(sp),cmap);
pcolor(...
    lon_OFAM, ...
    lat_OFAM, ...
    asal_now)
shading interp
hold on
caxis([asal_contours(1) asal_contours(end)]);
axis([110 152 -47 -31])
cbar = colorbar;
cbarrow
set(cbar, 'YTick',asal_contours(2:end-1));
Theta_contours = 10:1:21;
[c, h] = contour(...
    lon_OFAM, ...
    lat_OFAM, ...
    Theta_now, ...
    Theta_contours, ...
    'k', ...
    'linewidth',0.5);
clabel(c,h, ...
    'fontsize', font_size-1, ...
    'fontweight', 'bold')
[lon_mg, lat_mg]=meshgrid(lon_OFAM,lat_OFAM);
quiv_S = 5;
nn = 2;
u_now(u_now<0) = NaN;
quiver(...
    lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
    u_now(1:nn:end, 1:nn:end), ...
    v_now(1:nn:end, 1:nn:end), ...
    quiv_S, 'k');
% reference arrow goes here
title('OFAM surface')
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 11) 3rd subplot: CARS @ 400dbars
sp = 3;
axes(h_axes_sp(sp))
pres_now = 400;
pres_now_ind = pres == pres_now;
u_interp2_hr_now = ...
    interp2(lon_u_CARS, lat_u_CARS, ...
    u_CARS(:,:,pres_now_ind), lon_hr, lat_hr);
v_interp2_hr_now = ...
    interp2(lon_v_CARS, lat_v_CARS, ...
    v_CARS(:,:,pres_now_ind), lon_hr, lat_hr);
Theta_now = Theta_CARS(:,:,pres_now_ind);
asal_now = asal_CARS(:,:,pres_now_ind);
asal_contours = 34.78:0.02:35.16;
cmap_levels = length(asal_contours);
cmap = flipud(othercolor('RdBu11', cmap_levels-1));
colormap(h_axes_sp(sp),cmap);
pcolor(...
    lon_CARS, ...
    lat_CARS, ...
    asal_now)
shading interp
hold on
caxis([asal_contours(1) asal_contours(end)]);
axis([110 152 -47 -31])
cbar = colorbar;
set(cbar, 'visible', 'off');
Theta_contours = 9:0.2:12;
[c, h] = contour(...
    lon_CARS, ...
    lat_CARS, ...
    Theta_now, ...
    Theta_contours, ...
    'k', ...
    'linewidth',0.5);
clabel(c,h, ...
    'fontsize', font_size-1, ...
    'fontweight', 'bold')
[lon_mg, lat_mg]=meshgrid(lon_hr,lat_hr);
quiv_S = 5;
nn = 2;
u_interp2_hr_now(u_interp2_hr_now>0) = NaN;
quiver(...
    lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
    u_interp2_hr_now(1:nn:end, 1:nn:end), ...
    v_interp2_hr_now(1:nn:end, 1:nn:end), ...
    quiv_S, 'k');
% reference arrow goes here
title(['CARS ' num2str(pres_now) ' dbars'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 12) 4th subplot: OFAM @ 400fdb
sp = 4;
axes(h_axes_sp(sp))
pres_now = 400;
pres_now_ind = pres == pres_now;
u_now = u_OFAM(:,:,pres_now_ind);
v_now = v_OFAM(:,:,pres_now_ind);
Theta_now = Theta_OFAM(:,:,pres_now_ind);
asal_now = asal_OFAM(:,:,pres_now_ind);
asal_contours = 34.78:0.02:35.16;
cmap_levels = length(asal_contours);
cmap = flipud(othercolor('RdBu11', cmap_levels-1));
colormap(h_axes_sp(sp),cmap);
pcolor(...
    lon_OFAM, ...
    lat_OFAM, ...
    asal_now)
shading interp
hold on
caxis([asal_contours(1) asal_contours(end)]);
axis([110 152 -47 -31])
cbar = colorbar;
cbarrow
set(cbar, 'YTick',asal_contours(2:end-1));
Theta_contours = 9:0.2:12;
[c, h] = contour(...
    lon_OFAM, ...
    lat_OFAM, ...
    Theta_now, ...
    Theta_contours, ...
    'k', ...
    'linewidth',0.5);
clabel(c,h, ...
    'fontsize', font_size-1, ...
    'fontweight', 'bold')
[lon_mg, lat_mg]=meshgrid(lon_OFAM,lat_OFAM);
quiv_S = 5;
nn = 2;
u_now(u_now>0) = NaN;
quiver(...
    lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
    u_now(1:nn:end, 1:nn:end), ...
    v_now(1:nn:end, 1:nn:end), ...
    quiv_S, 'k');
% reference arrow goes here
title(['OFAM ' num2str(pres_now) ' dbars'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% print(gcf,'-r300','-dpng', ...
%     ['Dropbox/SACS_work/figures/f1_mean_TSV_maps/p_' ...
%     num2str(pres_now)])

export_fig(fig1, ['Dropbox/SACS_work/figures/f1_mean_TSV_maps/p_' ...
    num2str(pres_now)], ...
    '-m4', '-nocrop')
close

