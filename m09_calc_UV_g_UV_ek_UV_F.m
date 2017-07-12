%% calculate vertical integral of velocities
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_u_g'])
load([data_path 'SACS_data/aus8_v_g'])
load([data_path 'SACS_data/ERAInt'])
load([data_path 'SACS_data/aus8_rho'])


%%
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;
a = aus8_coor.a;
pi180 = aus8_coor.pi180;
depth = aus8_coor.depth;
depth_thkn = aus8_coor.depth_thkn;
depth_mid = aus8_coor.depth_mid;
Months = aus8_coor.Months;
U_mask = isnan(aus8_u_g.mean(:,:,1));
V_mask = isnan(aus8_v_g.mean(:,:,1));

aus8_coor.U_mask = U_mask;
aus8_coor.V_mask = V_mask;
save([data_path 'SACS_data/aus8_coor'], 'aus8_coor')


%% U_g
aus8_U_g.mean = NaN(length(lat_u), length(lon_u));
aus8_V_g.mean = NaN(length(lat_v), length(lon_v));
for t = 1 : 12
    aus8_U_g.(Months{t}) = ...
        NaN(length(lat_u), length(lon_u));
    aus8_V_g.(Months{t}) = ...
        NaN(length(lat_v), length(lon_v));
end

% integrate U along z. U_z is in m^2/s
u_bottom = NaN(length(lat_u), length(lon_u));
for m = 1 : length(lat_u)
    for n = 1 : length(lon_u)
        u_g_now = squeeze(aus8_u_g.mean(m,n,:));
        u_last_fini = find(isfinite(u_g_now), 1, 'last');
        u_first_nan = find(isnan(u_g_now), 1, 'first')-1;
        if isempty(u_last_fini)
            u_bottom(m,n) = 0;
        elseif isempty(u_first_nan)
            u_bottom(m,n) = -2000;
        else
            if u_last_fini == u_first_nan
                u_bottom(m,n) = depth(u_last_fini+1);
            else
                error('!')
            end
        end
        
        u_g_times_depth_thkn = u_g_now .* depth_thkn;
        aus8_U_g.mean(m,n) = nansum(u_g_times_depth_thkn);
        for t = 1 : 12
            u_g_now = squeeze(aus8_u_g.(Months{t})(m,n,:));
            u_g_times_depth_thkn = u_g_now .* depth_thkn;
            aus8_U_g.(Months{t})(m,n) = nansum(u_g_times_depth_thkn);
        end
    end
end
aus8_U_g.mean(U_mask) = NaN;
for t = 1 : 12
    aus8_U_g.(Months{t})(U_mask) = NaN;
end
aus8_coor.u_bottom = u_bottom;


% V_g
% integrate V along z. V_z is in m^2/s
v_bottom = NaN(length(lat_v), length(lon_v));
for m = 1 : length(lat_v)
    for n = 1 : length(lon_v)
        v_g_now = squeeze(aus8_v_g.mean(m,n,:));
        v_last_fini = find(isfinite(v_g_now), 1, 'last');
        v_first_nan = find(isnan(v_g_now), 1, 'first')-1;
        if isempty(v_last_fini)
            v_bottom(m,n) = 0;
        elseif isempty(v_first_nan)
            v_bottom(m,n) = -2000;
        else
            if v_last_fini == v_first_nan
                v_bottom(m,n) = depth(v_last_fini+1);
            else
                error('!')
            end
        end
        
        v_g_times_depth_thkn = v_g_now .* depth_thkn;
        aus8_V_g.mean(m,n) = nansum(v_g_times_depth_thkn);
        for t = 1 : 12
            v_g_now = squeeze(aus8_v_g.(Months{t})(m,n,:));
            v_g_times_depth_thkn = v_g_now .* depth_thkn;
            aus8_V_g.(Months{t})(m,n) = nansum(v_g_times_depth_thkn);
        end
    end
end
aus8_V_g.mean(V_mask) = NaN;
for t = 1 : 12
    aus8_V_g.(Months{t})(V_mask) = NaN;
end
aus8_coor.v_bottom = v_bottom;


%% U_ek V_ek
% interp Tau into corresponding uv grid
Tau_x_on_v.mean = interp2(lon_u, lat_v, ERAInt.Tau_x.mean, lon_v, lat_v);
Tau_y_on_u.mean = interp2(lon_u, lat_v, ERAInt.Tau_y.mean, lon_u, lat_u);
for t = 1 : 12
    Tau_x_on_v.(Months{t}) = ...
        interp2(lon_u, lat_v, ERAInt.Tau_x.(Months{t}), lon_v, lat_v);
    Tau_y_on_u.(Months{t}) = ...
        interp2(lon_u, lat_v, ERAInt.Tau_y.(Months{t}), lon_u, lat_u);
end

f_x = gsw_f(lat_u);
f_x_repmat = repmat(f_x, 1, length(lon_u));
f_y = gsw_f(lat_v);
f_y_repmat = repmat(f_y, 1, length(lon_v));

%
ekman_depth = -30;
ekman_layer_ind = depth >= ekman_depth;

%
rho_0 = nanmean(nanmean(nanmean(aus8_rho.mean(:,:,ekman_layer_ind))));

% Southern Hemisphere !
aus8_U_ek.mean = Tau_y_on_u.mean./f_x_repmat/rho_0;
aus8_V_ek.mean = -Tau_x_on_v.mean./f_y_repmat/rho_0;
aus8_U_ek.mean(U_mask) = NaN;
aus8_V_ek.mean(V_mask) = NaN;
for t = 1 : 12
    aus8_U_ek.(Months{t}) = Tau_y_on_u.(Months{t})./f_x_repmat/rho_0;
    aus8_V_ek.(Months{t}) = -Tau_x_on_v.(Months{t})./f_y_repmat/rho_0;
    aus8_U_ek.(Months{t})(U_mask) = NaN;
    aus8_V_ek.(Months{t})(V_mask) = NaN;
end


%% UV
aus8_U.mean = aus8_U_g.mean + aus8_U_ek.mean;
aus8_V.mean = aus8_V_g.mean + aus8_V_ek.mean;
aus8_U.mean(U_mask) = 0;
aus8_V.mean(V_mask) = 0;
for t = 1 : 12
    aus8_U.(Months{t}) = ...
        aus8_U_g.(Months{t}) + aus8_U_ek.(Months{t});
    aus8_V.(Months{t}) = ...
        aus8_V_g.(Months{t}) + aus8_V_ek.(Months{t});
    aus8_U.(Months{t})(U_mask) = 0;
    aus8_V.(Months{t})(V_mask) = 0;
end


%% F
dU.mean = (aus8_U.mean(:,2:end) - aus8_U.mean(:,1:end-1));
for t = 1 : 12
    dU.(Months{t}) = (aus8_U.(Months{t})(:,2:end) - ...
        aus8_U.(Months{t})(:,1:end-1));
end

dx = NaN(size(dU.mean));
for m = 1 : length(lat_u)
    dx_now = a * cos(lat_u(m) * pi180) .* ...
        (lon_u(2:end) - lon_u(1:end-1)) * pi180;
    
    dx(m,:) = dx_now;
end

lat_repmat = repmat(lat_v,1,length(lon_v));

dV.mean = ...
    aus8_V.mean(1:end-1,:) .* ...
    cos(lat_repmat(1:end-1,:) * pi180) - ...
    aus8_V.mean(2:end,:) .* ...
    cos(lat_repmat(2:end,:) * pi180);
for t = 1 : 12
    dV.(Months{t}) = ...
        aus8_V.(Months{t})(1:end-1,:) .* ...
        cos(lat_repmat(1:end-1,:) * pi180) - ...
        aus8_V.(Months{t})(2:end,:) .* ...
        cos(lat_repmat(2:end,:) * pi180);
end

dy = NaN(size(dV.mean));
for jj = 1 : length(lon_v)
    dy_now = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
    
    dy(:,jj) = dy_now;
end
lat_u_repmat = repmat(lat_u,1,length(lon_v));

aus8_F.mean = dU.mean./dx + 1./cos(lat_u_repmat * pi180).*dV.mean./dy;
for t = 1 : 12
    aus8_F.(Months{t}) = dU.(Months{t})./dx + ...
        1./cos(lat_u_repmat * pi180).*dV.(Months{t})./dy;
end

aus8_coor.F_mask = aus8_F.mean == 0;
save([data_path 'SACS_data/aus8_coor'], 'aus8_coor')


%%
save([data_path 'SACS_data/aus8_U_g_itgr'], 'aus8_U_g')
save([data_path 'SACS_data/aus8_V_g_itgr'], 'aus8_V_g')
save([data_path 'SACS_data/aus8_U_ek'], 'aus8_U_ek')
save([data_path 'SACS_data/aus8_V_ek'], 'aus8_V_ek')
save([data_path 'SACS_data/aus8_U'], 'aus8_U')
save([data_path 'SACS_data/aus8_V'], 'aus8_V')
save([data_path 'SACS_data/aus8_F'], 'aus8_F')

disp(['aus8_U_g_itgr aus8_V_g_itgr aus8_U_ek aus8_V_ek ' ...
    'aus8_U aus8_V aus8_F DONE'])

