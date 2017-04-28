%% calculate horizontal divergence
clearvars('-except', '*_path')

load aus8_ZD_method
a = aus8_ZD_method.a;
pi180 = aus8_ZD_method.pi180;
U = aus8_ZD_method.U;
lat_u = aus8_ZD_method.lat_u;
lon_u = aus8_ZD_method.lon_u;
V = aus8_ZD_method.V;
lat_v = aus8_ZD_method.lat_v;
lon_v = aus8_ZD_method.lon_v;
u_g_bottom_depth = aus8_ZD_method.u_g_bottom_depth;
v_g_bottom_depth = aus8_ZD_method.v_g_bottom_depth;


%% calculate horizontal divergence 
U(isnan(U)) = 0;
V(isnan(V)) = 0;

du = (U(:,2:end) - U(:,1:end-1));

dx = NaN(size(du));
for ii = 1 : length(lat_u)
    dx_now = a * cos(lat_u(ii) * pi180) .* ...
        (lon_u(2:end) - lon_u(1:end-1)) * pi180;
    
    dx(ii,:) = dx_now;
end

lat_repmat = repmat(lat_v,1,length(lon_v));
dv = V(1:end-1,:).*cos(lat_repmat(1:end-1,:) * pi180) - ...
    V(2:end,:).*cos(lat_repmat(2:end,:) * pi180);

dy = NaN(size(dv));
for jj = 1 : length(lon_v)
    dy_now = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
    
    dy(:,jj) = dy_now;
end

lat_u_repmat = repmat(lat_u,1,length(lon_v));
F = du./dx + 1./cos(lat_u_repmat * pi180).*dv./dy;

% Also get F's bottom depth for later
F_bottom_depth = NaN(size(F));
for ii = 1 : length(lat_u)
    for jj = 1 : length(lon_v)
        F_bottoms = ...
            [u_g_bottom_depth(ii,jj), u_g_bottom_depth(ii,jj+1), ...
            v_g_bottom_depth(ii,jj), v_g_bottom_depth(ii+1,jj)];
        max_depth_ind = find(F_bottoms == max(F_bottoms), 1, 'first');
        if isempty(max_depth_ind)
            continue
        else
            F_bottom_depth(ii,jj) = F_bottoms(max_depth_ind);
        end
    end
end

aus8_ZD_method.lat_F = lat_u;
aus8_ZD_method.lon_F = lon_v;
aus8_ZD_method.F = F;
aus8_ZD_method.F_bottom_depth = F_bottom_depth;


%%
save([cars_out_path 'aus8_ZD_method'], 'aus8_ZD_method')
disp(['aus8_ZD_method saved in ' ...
    cars_out_path 'aus8_ZD_method'])

