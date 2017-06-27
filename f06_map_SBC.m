%% fig 1: map of SBC currents
clearvars('-except', '*_path')

load aus8_ZD_method
load aus8_currents
lat = aus8_ZD_method.lat_phi;
lon = aus8_ZD_method.lon_phi;
pres = aus8_ZD_method.pres_mid;
depth_thicknesses = aus8_ZD_method.depth_thicknesses;
lat_u = aus8_ZD_method.lat_u;
lon_u = aus8_ZD_method.lon_u;
u_g_prime = aus8_ZD_method.u_g_prime;
lat_v = aus8_ZD_method.lat_v;
lon_v = aus8_ZD_method.lon_v;
v_g_prime = aus8_ZD_method.v_g_prime;
SBC_u_g_prime_ind = aus8_currents.SC.u_g_prime_ind;
SBC_v_g_prime_ind = aus8_currents.SC.v_g_prime_ind;
U_ek = aus8_ZD_method.U_ek;
V_ek = aus8_ZD_method.V_ek;
SC_lat_v_north = aus8_currents.SC.lat_v_north;
SC_lat_v_south = aus8_currents.SC.lat_v_south;
SC_lon_u_repelem = aus8_currents.SC.lon_u_repelem;
SC_lat_v_north_repelem = aus8_currents.SC.lat_v_north_repelem;
SC_lat_v_south_repelem = aus8_currents.SC.lat_v_south_repelem;


%%
u_dz = repmat(permute(depth_thicknesses, [3 2 1]), [size(U_ek), 1]);
u_g_prime_times_dz = u_g_prime .* u_dz;
SBC_U_g_prime = NaN(length(lat_u), length(lon_u));
for ii = 1 : length(lat_u)
    for jj = 1 : length(lon_u)
        u_g_prime_times_dz_now = squeeze(u_g_prime_times_dz(ii,jj,:));
        SBC_u_g_prime_ind_now = squeeze(SBC_u_g_prime_ind(ii,jj,:));
        SBC_U_g_prime(ii,jj) = ...
            nansum(u_g_prime_times_dz_now(SBC_u_g_prime_ind_now));
    end
end
SBC_U_g_prime(isnan(U_ek)) = NaN;

SBC_U_g_prime_mask = SBC_U_g_prime==0;
SBC_U_prime = SBC_U_g_prime + U_ek;
SBC_U_prime(SBC_U_g_prime==0) = 0;

v_dz = repmat(permute(depth_thicknesses, [3 2 1]), [size(V_ek), 1]);
v_g_prime_times_dz = v_g_prime .* v_dz;
SBC_V_g_prime = NaN(length(lat_v), length(lon_v));
for ii = 1 : length(lat_v)
    for jj = 1 : length(lon_v)
        v_g_prime_times_dz_now = v_g_prime_times_dz(ii,jj,:);
        SBC_v_g_prime_ind_now = squeeze(SBC_v_g_prime_ind(ii,jj,:));
        SBC_V_g_prime(ii,jj) = ...
            nansum(v_g_prime_times_dz(ii,jj,SBC_v_g_prime_ind_now));
    end
end
SBC_V_g_prime(SBC_V_g_prime==0) = NaN;
SBC_V_prime = SBC_V_g_prime + V_ek;

lon_hr = lon(1) : 0.1 : lon(end);
lat_hr = (lat(1) : -0.1 : lat(end))';

SBC_U_prime_interp2 = interp2(lon_u, lat_u, SBC_U_prime, ...
    lon, lat);
SBC_V_prime_interp2 = interp2(lon_v, lat_v, SBC_V_prime, ...
    lon, lat);


%% plot map
% 1) figure set-up
close all
font_size = 8;
fig1 = figure(1);
set(gcf, 'units', 'normalized', 'outerposition', [0 0 0.6 0.475]);
rowN = 1; colN = 1;
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
colormap(flipud(othercolor('Reds9')))

% 2) all data set-up
sp = 1;
axes(h_axes_sp(sp))
SBC_U_prime_now = SBC_U_prime;
SBC_U_prime_now(SBC_U_prime_now<0) = 0;

% 3) asal pcolor set-up
% new colormap routine !
levels = 11;
Reds = othercolor('Reds9', levels);
Reds(1,:) = [1 1 1];
magnif = 100;
Reds_cont = ...
    [0 0.05 0.1 0.5 1 5 10 20 50 100 200 300]*magnif;
Reds_cont_length = length(Reds_cont);
cmap_custom = cmapcust(Reds,Reds_cont);
colormap(cmap_custom);

% 4) plot asal pcolor
pcolor(...
    lon_u, ...
    lat_u, ...
    SBC_U_prime_now*magnif)
shading interp
hold on
caxis([Reds_cont(1) Reds_cont(end)]);
axis([114 148 -45 -33])
freezeColors

cmap = colormap(Reds);
Reds_linspace = ...
    linspace(Reds_cont(1), Reds_cont(end), Reds_cont_length);
cbar = colorbar;
set(cbar, 'YTick',Reds_linspace, 'YTickLabel',Reds_cont/magnif);

% 5) plot SBC boundaries
hold on
plot([SC_lon_u_repelem(1,1) SC_lon_u_repelem(1,1)], ...
    [SC_lat_v_north(1,1) SC_lat_v_south(1,1)], ...
    'g','linewidth',0.5,'linestyle','-')
plot(SC_lon_u_repelem(1,:), ...
    SC_lat_v_north_repelem(1,:), ...
    'g','linewidth',0.5,'linestyle','-')
plot(SC_lon_u_repelem(1,:), ...
    SC_lat_v_south_repelem(1,:), ...
    'g','linewidth',0.5,'linestyle','-')
plot([SC_lon_u_repelem(1,end) SC_lon_u_repelem(1,end)], ...
    [SC_lat_v_north(1,end) SC_lat_v_south(1,end)], ...
    'g','linewidth',0.5,'linestyle','-')

% 7) UV quiver set-up
[lon_mg, lat_mg]=meshgrid(lon,lat);
quiv_S = 5;
nn = 1;

% 8) plot directional U and V
quiver(...
    lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
    SBC_U_prime_interp2(1:nn:end, 1:nn:end), ...
    SBC_V_prime_interp2(1:nn:end, 1:nn:end), ...
    quiv_S, 'k');
% reference arrow goes here

% 9) title, grid, background and fonts
title('SBC current')
grid
set(gca,'xtick',115:2:147)
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
export_fig(fig1, [figures_path mfilename '/' scriptname(1:3) '_'], ...
    '-m4')
close
