%% calculate vertical integral of velocities
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_u_g'])
load([data_path 'SACS_data/aus8_v_g'])


%%
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
u_g = aus8_u_g.mean;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;
v_g = aus8_v_g.mean;
a = aus8_coor.a;
pi180 = aus8_coor.pi180;
depth = aus8_coor.depth;
depth_thkn = depth(2:end) - depth(1:end-1);
depth_mid = depth(1:end-1) + depth_thkn/2;


%% U
U_g = NaN(length(lat_u), length(lon_u));
% integrate U along z. U_z is in m^2/s
for ii = 1 : length(lat_u)
    for jj = 1 : length(lon_u)
        u_mid_now = squeeze(u_g_mid(ii,jj,:));
        
        u_mid_now_times_depth_thicknesses = ...
            u_mid_now .* depth_thicknesses;
        
        U_g(ii,jj) = nansum(u_mid_now_times_depth_thicknesses);
        
    end
end
u_mid_mask = isnan(u_g_mid(:,:,1));
U_g(u_mid_mask) = NaN;


% V
lat_v_length = length(lat_v);
lon_v_length = length(lon_v);
v_g_mid = NaN(lat_v_length, lon_v_length, pres_mid_length);
% get v at vertical halfway points
for ii = 1 : lat_v_length
    for jj = 1 : lon_v_length
        v_mid_now = interp1(pres,squeeze(v_g(ii,jj,:)),pres_mid);
        dx = a * cos(lat_v(ii) * pi180) .* ...
            (lon_u(jj+1) - lon_u(jj)) * pi180;
        
        v_g_mid(ii,jj,:) = interp1(pres,squeeze(v_g(ii,jj,:)),pres_mid);
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


% make topography mask
bottom_pres = 10 : 10 : 2000;
bottom_depth = gsw_z_from_p(-bottom_pres,lat_f_plane_p_to_depth)';

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


aus8_ZD_method.a = a;
aus8_ZD_method.pi180 = pi180;
aus8_ZD_method.lat_f_plane_p_to_depth = lat_f_plane_p_to_depth;

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


%%
save([cars_out_path 'aus8_ZD_method'], 'aus8_ZD_method')
disp(['aus8_ZD_method saved in ' ...
    cars_out_path 'aus8_ZD_method'])

