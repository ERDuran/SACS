%% calculate vertical integral of velocities
clearvars('-except', '*_path')

load aus8_geostrophy
lat_u = aus8_geostrophy.u.lat_u;
lon_u = aus8_geostrophy.u.lon_u;
u = aus8_geostrophy.u.u_0_HH_2000;
lat_v = aus8_geostrophy.v.lat_v;
lon_v = aus8_geostrophy.v.lon_v;
v = aus8_geostrophy.v.v_0_HH_2000;
a = aus8_geostrophy.a.value;
pi180 = aus8_geostrophy.a.pi180;
pres = aus8_geostrophy.pres;
pres_mid = (5 : 10 : 1995)';
pres_mid_length = length(pres_mid);


%% U
% depth on f-plane
lat_f_plane_p_to_depth = 40;
depth = gsw_z_from_p(-pres, lat_f_plane_p_to_depth);
depth_mid = gsw_z_from_p(-pres_mid, lat_f_plane_p_to_depth);
depth_thicknesses = depth(2:end) - depth(1:end-1);

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


% V
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

