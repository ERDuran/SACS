%% check that div(U_zd,V_zd) = 0
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

load aus8_ZD_method


%%
depth_thicknesses = aus8_ZD_method.depth_thicknesses;

u_g_prime = aus8_ZD_method.u_g_prime;

v_g_prime = aus8_ZD_method.v_g_prime;

a = aus8_ZD_method.a.value;
pi180 = aus8_ZD_method.a.pi180;
n_it = aus8_ZD_method.n_iteration;

lat_u = aus8_ZD_method.lat_u;
lon_u = aus8_ZD_method.lon_u;
lat_v = aus8_ZD_method.lat_v;
lon_v = aus8_ZD_method.lon_v;

U_ek = aus8_ZD_method.U_ek;
V_ek = aus8_ZD_method.V_ek;

F_nan = aus8_ZD_method.F == 0;


%% calculate vertical integration
lat_u_length = length(lat_u);
lon_u_length = length(lon_u);

U_g_prime = NaN(lat_u_length, lon_u_length);
% integrate U along z. U_z is in m^2/s
for ii = 1 : lat_u_length
    for jj = 1 : lon_u_length
        u_mid_now = squeeze(u_g_prime(ii,jj,:));
        
        u_mid_now_times_depth_thicknesses = ...
            u_mid_now .* depth_thicknesses;
        
        U_g_prime(ii,jj) = nansum(u_mid_now_times_depth_thicknesses);
        
    end
end

u_mid_mask = isnan(u_g_prime(:,:,1));
U_g_prime(u_mid_mask) = NaN;

lat_v_length = length(lat_v);
lon_v_length = length(lon_v);

V_g_prime = NaN(lat_v_length, lon_v_length);
% integrate V along z. V_z is in m^2/s
for ii = 1 : lat_v_length
    for jj = 1 : lon_v_length
        v_mid_now = squeeze(v_g_prime(ii,jj,:));
        
        v_mid_now_times_depth_thicknesses = ...
            v_mid_now .* depth_thicknesses;
        
        V_g_prime(ii,jj) = nansum(v_mid_now_times_depth_thicknesses);
        
    end
end

v_mid_mask = isnan(v_g_prime(:,:,1));
V_g_prime(v_mid_mask) = NaN;


%% Re-add the Ekman drift
U_prime = U_g_prime + U_ek;
V_prime = V_g_prime + V_ek;

U_prime(isnan(U_g_prime)) = 0;
V_prime(isnan(V_g_prime)) = 0;


%% calculate divergence
du = (U_prime(:,2:end) - U_prime(:,1:end-1));

dx = NaN(size(du));
for ii = 1 : length(lat_u)
    dx_now = a * cos(lat_u(ii) * pi180) .* ...
        (lon_u(2:end) - lon_u(1:end-1)) * pi180;
    
    dx(ii,:) = dx_now;
end

lat_v_repmat = repmat(lat_v,1,length(lon_v));
dv = V_prime(1:end-1,:).*cos(lat_v_repmat(1:end-1,:) * pi180) - ...
    V_prime(2:end,:).*cos(lat_v_repmat(2:end,:) * pi180);

dy = NaN(size(dv));
for jj = 1 : length(lon_v)
    dy_now = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
    
    dy(:,jj) = dy_now;
end

lat_u_repmat = repmat(lat_u,1,length(lon_v));
div_UV_prime = du./dx + 1./cos(lat_u_repmat * pi180).*dv./dy;
U_prime(isnan(U_g_prime)) = NaN;
V_prime(isnan(V_g_prime)) = NaN;
div_UV_prime(F_nan) = NaN;


%% 
close all
p15_check_uvzd


export_fig(['figures/p15_UV_d_prime_div/UV_d_prime_div_' num2str(n_it)], ...
    '-m3', '-nocrop');


%%
aus8_ZD_method.U_prime = U_prime;
aus8_ZD_method.V_prime = V_prime;
aus8_ZD_method.div_UV_prime = div_UV_prime;

save('cars_out/aus8_ZD_method', 'aus8_ZD_method')
disp('aus8_ZD_method saved in cars_out/aus8_ZD_method !')

