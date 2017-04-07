%% calculate vertical integral of velocities
% clc
path(pathdef)

% set up main directory
cd /home/z5100174/Desktop/MATLAB

% Add path to the data to be loaded
addpath(genpath('functions'))

% Add path to the data to be loaded
addpath cars_out
addpath programs/mean

% Add path to the gsw TEOS-10 library (including subfolders)
addpath(genpath('teos10_library'))

clear 
% vn is a function that can store the name of a variable
% hence vn(x) returns 'x' as a string.
vn = @(x) inputname(1);

load aus8_geostrophy


%%
lat_u = aus8_geostrophy.u.lat_u;
lon_u = aus8_geostrophy.u.lon_u;
u = aus8_geostrophy.u.u_0_HH;

lat_v = aus8_geostrophy.v.lat_v;
lon_v = aus8_geostrophy.v.lon_v;
v = aus8_geostrophy.v.v_0_HH;

cst_lat = aus8_geostrophy.f.cst_lat;
a = aus8_geostrophy.a.value;
pi180 = aus8_geostrophy.a.pi180;

%
pres = aus8_geostrophy.pres;
pres_mid = (5 : 10 : 1995)';
pres_mid_length = length(pres_mid);

%
depth = gsw_z_from_p(-pres, cst_lat);
depth_mid = gsw_z_from_p(-pres_mid, cst_lat);
depth_thicknesses = depth(2:end) - depth(1:end-1);

%
% smooth_weigth = aus8_geostrophy.smooth_weigth;


%%
lat_u_length = length(lat_u);
lon_u_length = length(lon_u);


u_g_mid = NaN(lat_u_length, lon_u_length, pres_mid_length);
% get u at vertical halfway points
for ii = 1 : lat_u_length
    for jj = 1 : lon_u_length
        u_mid_now = interp1(pres,squeeze(u(ii,jj,:)),pres_mid);
        dy = a * (lat_v(ii) - lat_v(ii+1)) * pi180;
        
        u_g_mid(ii,jj,:) = u_mid_now; 
    end
end


U_g = NaN(lat_u_length, lon_u_length);
% integrate U along z. U_z is in m^2/s
for ii = 1 : lat_u_length
    for jj = 1 : lon_u_length
        u_mid_now = squeeze(u_g_mid(ii,jj,:));
        
        u_mid_now_times_depth_thicknesses = ...
            u_mid_now .* depth_thicknesses;
        
        U_g(ii,jj) = nansum(u_mid_now_times_depth_thicknesses);
        
    end
end

u_mid_mask = isnan(u_g_mid(:,:,1));
U_g(u_mid_mask) = NaN;


%%
lat_v_length = length(lat_v);
lon_v_length = length(lon_v);


v_g_mid = NaN(lat_v_length, lon_v_length, pres_mid_length);
% get v at vertical halfway points
for ii = 1 : lat_v_length
    for jj = 1 : lon_v_length
        v_mid_now = interp1(pres,squeeze(v(ii,jj,:)),pres_mid);
        dx = a * cos(lat_v(ii) * pi180) .* ...
            (lon_u(jj+1) - lon_u(jj)) * pi180;
        
        v_g_mid(ii,jj,:) = interp1(pres,squeeze(v(ii,jj,:)),pres_mid);
    end
end


V_g = NaN(lat_v_length, lon_v_length);
% integrate V along z. V_z is in m^2/s
for ii = 1 : lat_v_length
    for jj = 1 : lon_v_length
        v_mid_now = squeeze(v_g_mid(ii,jj,:));
        
        v_mid_now_times_depth_thicknesses = ...
            v_mid_now .* depth_thicknesses;
        
        V_g(ii,jj) = nansum(v_mid_now_times_depth_thicknesses);
        
    end
end

v_mid_mask = isnan(v_g_mid(:,:,1));
V_g(v_mid_mask) = NaN;


%% make topography mask
bottom_pres = 10 : 10 : 2000;
bottom_depth = gsw_z_from_p(-bottom_pres,cst_lat)';

u_bottom_depth = NaN(length(lat_u), length(lon_u));
for ii = 1 : length(lat_u)
    for jj = 1 : length(lon_u)
        bottom_ind = find(isfinite(u_g_mid(ii,jj,:)), 1, 'last');
        
        if isempty(bottom_ind)
            continue
        end
        
        u_bottom_depth(ii,jj) = bottom_depth(bottom_ind);
    end
end


v_bottom_depth = NaN(length(lat_v), length(lon_v));
for ii = 1 : length(lat_v)
    for jj = 1 : length(lon_v)
        bottom_ind = find(isfinite(v_g_mid(ii,jj,:)), 1, 'last');
        
        if isempty(bottom_ind)
            continue
        end
        
        v_bottom_depth(ii,jj) = bottom_depth(bottom_ind);
    end
end


%% svdkjbj
u_plot = u_g_mid;
v_plot = v_g_mid;

pres_mid_1 = 5;
pres_mid_2 = 255;

close all
p10_p16_plot_cross_sections

export_fig(fig1, ['figures/p10_plot_ugvg/merid_u5'], ...
    '-m3', '-nocrop')

% export_fig(fig2, ['figures/p10_plot_ugvg/merid_v'], ...
%     '-m3', '-nocrop')

export_fig(fig3, ['figures/p10_plot_ugvg/merid_maps5'], ...
    '-m3', '-nocrop')
% close all


%%
aus8_ZD_method.a.value = a;
aus8_ZD_method.f.cst_lat = cst_lat;
aus8_ZD_method.a.pi180 = pi180;

aus8_ZD_method.pres_mid = pres_mid;

aus8_ZD_method.depth = depth;
aus8_ZD_method.depth_mid = depth_mid;
aus8_ZD_method.depth_thicknesses = depth_thicknesses;

aus8_ZD_method.lat_u = lat_u;
aus8_ZD_method.lon_u = lon_u;
aus8_ZD_method.u_g = u_g_mid;
aus8_ZD_method.u_g_bottom_depth = u_bottom_depth;
aus8_ZD_method.U_g = U_g;

aus8_ZD_method.lat_v = lat_v;
aus8_ZD_method.lon_v = lon_v;
aus8_ZD_method.v_g = v_g_mid;
aus8_ZD_method.v_g_bottom_depth = v_bottom_depth;
aus8_ZD_method.V_g = V_g;

save('cars_out/aus8_ZD_method', 'aus8_ZD_method')
disp('aus8_ZD_method saved in cars_out/aus8_ZD_method !')

