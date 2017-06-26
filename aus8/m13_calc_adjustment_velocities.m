%% calculate adjustment velocities
clearvars('-except', '*_path')

load aus8_ZD_method
lat_f_plane_p_to_depth = aus8_ZD_method.lat_f_plane_p_to_depth;
pres_mid = aus8_ZD_method.pres_mid;
depth_thicknesses = aus8_ZD_method.depth_thicknesses;
depth_mid = aus8_ZD_method.depth_mid;
lat_u = aus8_ZD_method.lat_u;
lon_u = aus8_ZD_method.lon_u;
u_g = aus8_ZD_method.u_g;
lat_v = aus8_ZD_method.lat_v;
lon_v = aus8_ZD_method.lon_v;
v_g = aus8_ZD_method.v_g;
U_d = aus8_ZD_method.U_d;
V_d = aus8_ZD_method.V_d;


%% make topography mask
bottom_pres = 10 : 10 : 2000;
bottom_depth = gsw_z_from_p(-bottom_pres,lat_f_plane_p_to_depth)';

u_bottom_depth = NaN(size(U_d));
for ii = 1 : length(lat_u)
    for jj = 1 : length(lon_u)
        bottom_ind = find(isfinite(u_g(ii,jj,:)), 1, 'last');
        
        if isempty(bottom_ind)
            continue
        end
        
        u_bottom_depth(ii,jj) = bottom_depth(bottom_ind);
    end
end


v_bottom_depth = NaN(size(V_d));
for ii = 1 : length(lat_v)
    for jj = 1 : length(lon_v)
        bottom_ind = find(isfinite(v_g(ii,jj,:)), 1, 'last');
        
        if isempty(bottom_ind)
            continue
        end
        
        v_bottom_depth(ii,jj) = bottom_depth(bottom_ind);
    end
end

%
u_d = U_d ./ u_bottom_depth;
u_d_repmat = repmat(u_d,1,1,length(pres_mid));
u_g_prime = u_g - u_d_repmat;

v_d = V_d ./ v_bottom_depth;
v_d_repmat = repmat(v_d,1,1,length(pres_mid));
v_g_prime = v_g - v_d_repmat;

aus8_ZD_method.u_g_prime = u_g_prime;
aus8_ZD_method.u_g_prime_bottom_depth = u_bottom_depth;
aus8_ZD_method.v_g_prime = v_g_prime;
aus8_ZD_method.v_g_prime_bottom_depth = v_bottom_depth;


%%
save([cars_out_path 'aus8_ZD_method'], 'aus8_ZD_method')
disp(['aus8_ZD_method saved in ' ...
    cars_out_path 'aus8_ZD_method'])

