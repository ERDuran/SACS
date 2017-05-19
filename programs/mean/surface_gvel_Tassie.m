%%
clearvars('-except', '*_path')

load aus8_ZD_method
load aus8_geostrophy

lat_u = aus8_ZD_method.lat_u;
lon_u = aus8_ZD_method.lon_u;
u_g = aus8_geostrophy.u.u_0_HH_2000;
lat_v = aus8_ZD_method.lat_v;
lon_v = aus8_ZD_method.lon_v;
v_g = aus8_geostrophy.v.v_0_HH_2000;
U_ek = aus8_ZD_method.U_ek;
V_ek = aus8_ZD_method.V_ek;


%%
U_prime = (u_g(:,:,1) + U_ek/5)*3.6;
V_prime = (v_g(:,:,1) + V_ek/40)*3.6;
% U_prime = (u_g_prime(:,:,1) + 0)*3.6;
% V_prime = (v_g_prime(:,:,1) + 0)*3.6;

U_interp2 = interp2(...
    lon_u, lat_u, U_prime, lon_v, lat_u);
V_interp2 = interp2(...
    lon_v, lat_v, V_prime, lon_v, lat_u);

speed_prime = sqrt(U_interp2.^2 + V_interp2.^2);
% speed_prime(isnan(speed_prime)) = 0;


%%
% 1) figure set-up
close all
font_size = 8;
fig1 = figure(1);
set(gcf, 'units', 'normalized', 'outerposition', [0 0 0.6 0.8]);
rowN = 1; colN = 1;
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = - 0.04; % gap width between subplots
gap_h = 0.04; % gap height between subplots
marg_b = 0.04; % bottom margin
marg_t = 0.04; % top margin
marg_l = 0.02; % left margin
marg_r = 0.02; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length
h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

% 2) all data set-up
sp = 1;
axes(h_axes_sp(sp))

% 3) asal pcolor set-up
% new colormap routine !
magnif = 100;
cmap2_cont = ...
    [0 0.05 0.1 0.2 0.3 0.4 0.5 0.75 1];
% cmap1_cont = -fliplr(cmap2_cont);
levels = length(cmap2_cont)-1;

% cmap1 = flipud(othercolor('Blues9', levels));
% cmap1(end,:) = [1 1 1];
% cmap2 = othercolor('Reds9', levels);
cmap2 = othercolor('Blues9', levels);
cmap2(1,:) = [1 1 1];
% cmaps = [cmap1; cmap2];

% cmaps_cont = [cmap1_cont cmap2_cont(2:end)];
% cmaps_cont_length = length(cmaps_cont);
cmap2_cont_length = length(cmap2_cont);
% cmaps_custom = cmapcust(cmaps,cmaps_cont);
cmaps_custom = cmapcust(cmap2,cmap2_cont*magnif);
colormap(h_axes_sp(sp), cmaps_custom);

% 4) plot asal pcolor
pcolor(...
    lon_v-1/16, ...
    lat_u+1/16, ...
    speed_prime)
shading interp
% freezeColors

hold on

caxis([cmap2_cont(1) cmap2_cont(end)]);
lon_min = 140; lon_max = 152; lat_min = -45; lat_max = -38;
axis([lon_min lon_max lat_min lat_max])
%
% cmap = colormap(cmap2);
% cmaps_linspace = ...
%     linspace(cmaps_cont(1), cmaps_cont(end), cmaps_cont_length);
% cmaps_linspace = ...
%     linspace(cmap2_cont(1), cmap2_cont(end), cmap2_cont_length);
cbar = colorbar;
% set(cbar, 'YTick',cmaps_linspace, 'YTickLabel',cmap2_cont);


% 7) UV quiver set-up
[lon_mg, lat_mg]=meshgrid(lon_v,lat_u);
qscale = 50/magnif  ; % scaling factor for all vectors
nn = 1;

% 8) plot directional U and V
h = quiver(...
    lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
    U_interp2(1:nn:end, 1:nn:end), ...
    V_interp2(1:nn:end, 1:nn:end), ...
    0, 'k');
hU = get(h,'UData') ;
hV = get(h,'VData') ;
set(h,'UData',qscale*hU,'VData',qscale*hV)
% reference vector at the end

% 9) title, grid, background and fonts
title(['Mean velocity at the surface in km/h'])
grid
set(gca,'xtick',lon_min:2:lon_max)
set(gca,'ytick',lat_min:1:lat_max)
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size)
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% reference vector
ref_vec = 0.5;
n_lon = 2;
lon_length = (lon_max - lon_min);
lon_x_ratio = x_sp / lon_length;
lon_x_ratio = lon_x_ratio*n_lon;
n_lat = 1;
lat_length = (lat_max - lat_min);
lat_y_ratio = y_sp / lat_length;
lat_y_ratio = lat_y_ratio*n_lat;
axes('Position', ...
    [marg_l+x_sp-lon_x_ratio-0.06 ...
    marg_b+y_sp-lat_y_ratio ...
    lon_x_ratio lat_y_ratio], ...
    'layer', 'top');
ref_x = 0:1/8:1*n_lon;
ref_y = 0:1/8:1*n_lat;
[ref_x_mg, ref_y_mg] = meshgrid(ref_x, ref_y);
ref_U = zeros(size(ref_x_mg));
ref_U(4*n_lat,4*n_lon) = ref_vec;
ref_V = zeros(size(ref_x_mg));

h = quiver(...
    ref_x, ref_y, ...
    ref_U, ...
    ref_V, ...
    0, 'k');
hU = get(h,'UData') ;
hV = get(h,'VData') ;
set(h,'UData',qscale*hU,'VData',qscale*hV)
set(gca, 'xtick', '')
set(gca, 'ytick', '')
text(ref_x(4*n_lon), ref_y(4*n_lat+2), [num2str(ref_vec) ' km/h'], ...
    'fontsize',font_size)


outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
export_fig(fig1, [figures_path mfilename '/' scriptname(1:3) '_'], ...
    '-m4')
close