%% fig 1: map at the surface of T in colours S in contours and V in arrows
clc
path(pathdef)

% set up main directory
cd ~
% cd /home/z5100174/Desktop/MATLAB

% Add path to the data to be loaded
addpath(genpath('Dropbox/SACS_work'))

clear 
% vn is a function that can store the name of a variable
% hence vn(x) returns 'x' as a string.
vn = @(x) inputname(1);

load aus8_TEOS10
load aus8_geostrophy


%% load stuff
lat = aus8_TEOS10.lat;
lon = aus8_TEOS10.lon;
pres = aus8_TEOS10.pres;
Theta = aus8_TEOS10.Theta.mean;
asal = aus8_TEOS10.asal.mean;

lat_u = aus8_geostrophy.u.lat_u;
lon_u = lon;
u = aus8_geostrophy.u.u_0_HH;

lat_v = lat;
lon_v = aus8_geostrophy.v.lon_v;
v = aus8_geostrophy.v.v_0_HH;


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
set(gcf, 'units', 'normalized', 'outerposition', [0 0 0.8 1]);
rowN = 2; colN = 2;
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (flipud(repmat((1:rowN)',1,colN)))';
gap_w = 0.02; % gap width between subplots
gap_h = 0.06; % gap height between subplots
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
axes(h_axes_sp(1))
pres_now = 0;
pres_now_ind = pres == pres_now;
lon_hr = lon(1) : 0.1 : lon(end);
lat_hr = (lat(1) : -0.1 : lat(end))';
u_interp2_hr_now = ...
    interp2(lon_u, lat_u, u(:,:,pres_now_ind), lon_hr, lat_hr);
v_interp2_hr_now = ...
    interp2(lon_v, lat_v, v(:,:,pres_now_ind), lon_hr, lat_hr);
Theta_now = Theta(:,:,pres_now_ind);
asal_now = asal(:,:,pres_now_ind);

% 3) asal pcolor set-up
asal_contours = 34.7:0.1:36.5;
cmap_levels = length(asal_contours);
cmap = flipud(othercolor('RdBu11', cmap_levels-1));
colormap(h_axes_sp(1), cmap);

% 4) plot asal pcolor
pcolor(...
    lon, ...
    lat, ...
    asal_now)
shading interp
hold on
caxis([asal_contours(1) asal_contours(end)]);
axis([110 152 -47 -31])
cbar = colorbar;
cbarrow
set(cbar, 'YTick',asal_contours(2:end-1));

% 5) Theta contours set-up
Theta_contours = 10:1:21;

% 6) plot Theta labelled contours
[c, h] = contour(...
    lon, ...
    lat, ...
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

% 10) 2nd subplot: OFAM @ surface


% 11) 3rd subplot: CARS @ 400dbars
axes(h_axes_sp(3))
pres_now = 400;
pres_now_ind = pres == pres_now;
u_interp2_hr_now = ...
    interp2(lon_u, lat_u, u(:,:,pres_now_ind), lon_hr, lat_hr);
v_interp2_hr_now = ...
    interp2(lon_v, lat_v, v(:,:,pres_now_ind), lon_hr, lat_hr);
Theta_now = Theta(:,:,pres_now_ind);
asal_now = asal(:,:,pres_now_ind);
asal_contours = 34.78:0.02:35.16;
cmap_levels = length(asal_contours);
cmap = flipud(othercolor('RdBu11', cmap_levels-1));
colormap(h_axes_sp(3),cmap);
pcolor(...
    lon, ...
    lat, ...
    asal_now)
shading interp
hold on
caxis([asal_contours(1) asal_contours(end)]);
axis([110 152 -47 -31])
cbar = colorbar;
cbarrow
set(cbar, 'YTick',asal_contours(2:end-1));
Theta_contours = 9:0.2:11;
[c, h] = contour(...
    lon, ...
    lat, ...
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

% 4th subplot: OFAM @ 400fdb


% print(gcf,'-r300','-dpng', ...
%     ['Dropbox/SACS_work/figures/f1_mean_TSV_maps/p_' ...
%     num2str(pres_now)])

export_fig(fig1, ['Dropbox/SACS_work/figures/f1_mean_TSV_maps/p_' ...
    num2str(pres_now)], ...
    '-m4', '-nocrop')
close





