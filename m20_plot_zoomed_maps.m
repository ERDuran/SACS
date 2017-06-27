%% define current boundaries
clc
path(pathdef)

% set up main directory
cd /home/z5100174/Desktop/MATLAB
% cd D:/

% Add path to the data to be loaded
addpath(genpath('functions'))

% Add path to the data to be loaded
addpath cars_out

% Add path to the gsw TEOS-10 library (including subfolders)
addpath(genpath('teos10_library'))

clear 
% vn is a function that can store the name of a variable
% hence vn(x) returns 'x' as a string.
vn = @(x) inputname(1);

load aus8_ZD_method
load aus8_streamfunction
load aus8_geostrophy


%% load data
% note since I want to plot a map using quiver,
% u and v must be on the same grid (T and S grid)
% therefore this will ruin the zero divergence property
% which is ok here because we just want to look at maps of
% geostrophic velocity fields

pres = (0:10:2000)';

pres_mid = aus8_ZD_method.pres_mid;

lat_u = aus8_ZD_method.lat_u;
lon_u = aus8_ZD_method.lon_u;
u_g_prime = aus8_ZD_method.u_g_prime;

lat_v = aus8_ZD_method.lat_v;
lon_v = aus8_ZD_method.lon_v;
v_g_prime = aus8_ZD_method.v_g_prime;

lat_dynh = aus8_streamfunction.lat;
lon_dynh = aus8_streamfunction.lon;
dynh = aus8_streamfunction.dynh.p_ref_0;

lat_F = aus8_ZD_method.lat_F;
lon_F = aus8_ZD_method.lon_F;
div_UV_prime = aus8_ZD_method.div_UV_prime;

F = aus8_ZD_method.F;
Lap_phi = aus8_ZD_method.Lap_phi;

%
dynh_mid = NaN(length(lat_dynh), length(lon_dynh), length(pres_mid));
for ii = 1 : length(lat_dynh)
    for jj = 1 : length(lon_dynh)
        dynh_mid(ii,jj,:) = ...
            interp1(pres, squeeze(dynh(ii,jj,:)), pres_mid);
        
    end
end


%% plots 
pres_mid_now = 5;
p_ind = find(pres_mid == pres_mid_now);

data.u_g = aus8_ZD_method.u_g(:,:,p_ind);
data.U_g = aus8_ZD_method.U_g;
data.U_ek = aus8_ZD_method.U_ek;
data.U = aus8_ZD_method.U;
data.U_d = aus8_ZD_method.U_d;
data.u_g_prime = u_g_prime(:,:,p_ind);
data.U_prime = aus8_ZD_method.U_prime;

data.v_g = aus8_ZD_method.v_g(:,:,p_ind);
data.V_g = aus8_ZD_method.V_g;
data.V_ek = aus8_ZD_method.V_ek;
data.V = aus8_ZD_method.V;
data.V_d = aus8_ZD_method.V_d;
data.v_g_prime = v_g_prime(:,:,p_ind);
data.V_prime = aus8_ZD_method.V_prime;

data.dynh = dynh_mid(:,:,p_ind);
data.div_UV_prime = div_UV_prime;
data.F = F;
data.Lap_phi = Lap_phi;

% define surface current boundaries
% longitude of western boundary of the control volume
ALLC_west_lon_u = 115;
% index of longitude of western boundary of the control volume
ALLC_west_lon_u_ind = find(lon_u == ALLC_west_lon_u);

% longitude of eastern boundary of the control volume
ALLC_east_lon_u = 146;
% index of longitude of eastern boundary of the control volume
ALLC_east_lon_u_ind = find(lon_u == ALLC_east_lon_u);

% length of all the u longitudes within the west and east boundaries
ALLC_lon_u_length = length(...
    ALLC_west_lon_u_ind : ALLC_east_lon_u_ind);

%
lon_u_ind = ...
    find(lon_u >= ALLC_west_lon_u & lon_u <= ALLC_east_lon_u);
lon_v_ind = ...
    find(lon_v >= ALLC_west_lon_u & lon_v <= ALLC_east_lon_u);
lon_dynh_ind_raw = ...
    find(lon_dynh >= ALLC_west_lon_u & lon_dynh <= ALLC_east_lon_u);
lon_dynh_ind = ...
    [lon_dynh_ind_raw(1)-1, lon_dynh_ind_raw, lon_dynh_ind_raw(end)+1];

%
lon_plot_west = 114;
lon_plot_east = 150;
lon_plot_ind =  ...
    find(lon_dynh >= lon_plot_west & lon_dynh <= lon_plot_east);

%
lat_plot_north = -31;
lat_plot_south = -46;

%
lat_dynh_ind = ...
    find(lat_dynh <= lat_plot_north & lat_dynh >= lat_plot_south);
lat_plot_ind = lat_dynh_ind;


% plot maps at some depths
close all
fig1 = figure(1);

set(gcf,'units','normalized','outerposition',[0 0 2 2])

rowN = 1; colN = 1;
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

axes(h_axes_sp(1))

% new colormap routine !
levels = 11;
Reds = othercolor('Reds9', levels);
Blues = othercolor('Blues9', levels);
[Reds(1,:), Blues(1,:)] = deal([1 1 1]);
Blues = flipud(Blues);
BluesReds = [Blues; Reds];
magnif = 1000;

dn = 'div_UV_prime';
dn_units = 'm^2/s';

if strcmp(dn_units, 'm^2/s')
    Reds_cont = [0 0.05 0.1 0.5 1 5 10 20 50 100 200 300]*magnif;
else
    Reds_cont = [0 0.005 0.01 0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5]*magnif;
end

Blues_cont = -fliplr(Reds_cont);

BluesReds_cont = [Blues_cont, Reds_cont(2:end)];
BluesReds_cont_length = length(BluesReds_cont);

cmap_custom = cmapcust(BluesReds,BluesReds_cont);
colormap(cmap_custom);

h_pcolor = pcolor(...
    lon_F(lon_plot_ind)-0.0625, ...
    lat_F(lat_plot_ind)+0.0625, ...
    data.(dn)(lat_plot_ind,lon_plot_ind)*magnif);

% shading flat
caxis([Blues_cont(1) Reds_cont(end)]);

freezeColors

hold on
nn = 1;
quiv_S = 10;
[lon_mg, lat_mg]=meshgrid(lon_u,lat_u);
v_g_0 = zeros(size(u_g_prime));
quiver(...
    lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
    u_g_prime(1:nn:end, 1:nn:end, p_ind), ...
    v_g_0(1:nn:end, 1:nn:end, p_ind), ...
    quiv_S, 'r');

hold on
nn = 1;
quiv_S = 10;
[lon_mg, lat_mg]=meshgrid(lon_v,lat_v);
u_g_0 = zeros(size(v_g_prime));
quiver(...
    lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
    u_g_0(1:nn:end, 1:nn:end, p_ind), ...
    v_g_prime(1:nn:end, 1:nn:end, p_ind), ...
    quiv_S, 'b');

hold on
[lat_u_g_prime_NaN_ind, lon_u_g_prime_NaN_ind] = ...
    find(isnan(u_g_prime(:,:,p_ind)));

lat_u_g_prime_NaN = lat_u(lat_u_g_prime_NaN_ind);
lon_u_g_prime_NaN = lon_u(lon_u_g_prime_NaN_ind);

scatter(lon_u_g_prime_NaN, lat_u_g_prime_NaN, 4, ...
    'o', 'k');

hold on
[lat_v_g_prime_NaN_ind, lon_v_g_prime_NaN_ind] = ...
    find(isnan(v_g_prime(:,:,p_ind)));

lat_v_g_prime_NaN = lat_v(lat_v_g_prime_NaN_ind);
lon_v_g_prime_NaN = lon_v(lon_v_g_prime_NaN_ind);

scatter(lon_v_g_prime_NaN, lat_v_g_prime_NaN, 4, ...
    'o', 'k');


if strcmp(dn_units, 'm^2/s')
    title([dn ' (' dn_units ')'])
else
    title([dn ' (' dn_units ')' ' at ' num2str(pres_mid_now) ' dbar'])
end

cmap = colormap(BluesReds);
BluesReds_linspace = ...
    linspace(Blues_cont(1), Reds_cont(end), BluesReds_cont_length);
cbar = colorbar;
% cbarrow
set(cbar, 'YTick',BluesReds_linspace, 'YTickLabel',BluesReds_cont/magnif);

% grid 
font_size = 10;
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')


export_fig(fig1, ['figures/m15_plot_zoomed_maps/' dn '_' ...
    num2str(pres_mid_now)], ...
    '-m4', '-nocrop')

close

