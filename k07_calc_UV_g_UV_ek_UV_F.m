%% calculate vertical integral of velocities
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/KDau_u_g'])
load([data_path 'SACS_data/KDau_v_g'])
load([data_path 'SACS_data/KDau_tau'])
load([data_path 'SACS_data/KDau_rho'])


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
U_mask_KDau_dum = isnan(KDau_u_g.mean(:,:,1));
V_mask_KDau_dum = isnan(KDau_v_g.mean(:,:,1));
tau_x_mask_KDau = isnan(KDau_tau.tau_x.mean);
tau_y_mask_KDau = isnan(KDau_tau.tau_y.mean);
U_mask_KDau = U_mask_KDau_dum | tau_y_mask_KDau;
V_mask_KDau = V_mask_KDau_dum | tau_x_mask_KDau;

aus8_coor.U_mask_KDau = U_mask_KDau;
aus8_coor.V_mask_KDau = V_mask_KDau;
save([data_path 'SACS_data/aus8_coor'], 'aus8_coor')


%% U_g
KDau_U_g.mean = NaN(length(lat_u), length(lon_u));
KDau_V_g.mean = NaN(length(lat_v), length(lon_v));
for t = 1 : 12
    KDau_U_g.(Months{t}) = ...
        NaN(length(lat_u), length(lon_u));
    KDau_V_g.(Months{t}) = ...
        NaN(length(lat_v), length(lon_v));
end

% integrate U along z. U_z is in m^2/s
u_bottom_KDau = NaN(length(lat_u), length(lon_u));
for m = 1 : length(lat_u)
    for n = 1 : length(lon_u)
        u_g_now = squeeze(KDau_u_g.mean(m,n,:));
        u_last_fini = find(isfinite(u_g_now), 1, 'last');
        u_first_nan = find(isnan(u_g_now), 1, 'first')-1;
        if isempty(u_last_fini)
            u_bottom_KDau(m,n) = 0;
        elseif isempty(u_first_nan)
            u_bottom_KDau(m,n) = -2000;
        else
            if u_last_fini == u_first_nan
                u_bottom_KDau(m,n) = depth(u_last_fini+1);
            else
                error('!')
            end
        end
        
        u_g_times_depth_thkn = u_g_now .* depth_thkn;
        KDau_U_g.mean(m,n) = nansum(u_g_times_depth_thkn);
        for t = 1 : 12
            u_g_now = squeeze(KDau_u_g.(Months{t})(m,n,:));
            u_g_times_depth_thkn = u_g_now .* depth_thkn;
            KDau_U_g.(Months{t})(m,n) = nansum(u_g_times_depth_thkn);
        end
    end
end
KDau_U_g.mean(U_mask_KDau) = NaN;
for t = 1 : 12
    KDau_U_g.(Months{t})(U_mask_KDau) = NaN;
end
aus8_coor.u_bottom_KDau = u_bottom_KDau;


% V_g
% integrate V along z. V_z is in m^2/s
v_bottom_KDau = NaN(length(lat_v), length(lon_v));
for m = 1 : length(lat_v)
    for n = 1 : length(lon_v)
        v_g_now = squeeze(KDau_v_g.mean(m,n,:));
        v_last_fini = find(isfinite(v_g_now), 1, 'last');
        v_first_nan = find(isnan(v_g_now), 1, 'first')-1;
        if isempty(v_last_fini)
            v_bottom_KDau(m,n) = 0;
        elseif isempty(v_first_nan)
            v_bottom_KDau(m,n) = -2000;
        else
            if v_last_fini == v_first_nan
                v_bottom_KDau(m,n) = depth(v_last_fini+1);
            else
                error('!')
            end
        end
        
        v_g_times_depth_thkn = v_g_now .* depth_thkn;
        KDau_V_g.mean(m,n) = nansum(v_g_times_depth_thkn);
        for t = 1 : 12
            v_g_now = squeeze(KDau_v_g.(Months{t})(m,n,:));
            v_g_times_depth_thkn = v_g_now .* depth_thkn;
            KDau_V_g.(Months{t})(m,n) = nansum(v_g_times_depth_thkn);
        end
    end
end
KDau_V_g.mean(V_mask_KDau) = NaN;
for t = 1 : 12
    KDau_V_g.(Months{t})(V_mask_KDau) = NaN;
end
aus8_coor.v_bottom_KDau = v_bottom_KDau;


%% U_ek V_ek
f_x = gsw_f(lat_u);
f_x_repmat = repmat(f_x, 1, length(lon_u));
f_y = gsw_f(lat_v);
f_y_repmat = repmat(f_y, 1, length(lon_v));

%
ekman_depth = -30;
ekman_layer_ind = depth >= ekman_depth;

%
rho_0 = nanmean(nanmean(nanmean(KDau_rho.mean(:,:,ekman_layer_ind))));

% Southern Hemisphere !
KDau_U_ek.mean = KDau_tau.tau_y.mean./f_x_repmat/rho_0;
KDau_V_ek.mean = -KDau_tau.tau_x.mean./f_y_repmat/rho_0;
KDau_U_ek.mean(U_mask_KDau) = NaN;
KDau_V_ek.mean(V_mask_KDau) = NaN;
for t = 1 : 12
    KDau_U_ek.(Months{t}) = KDau_tau.tau_y.(Months{t})./f_x_repmat/rho_0;
    KDau_V_ek.(Months{t}) = -KDau_tau.tau_x.(Months{t})./f_y_repmat/rho_0;
    KDau_U_ek.(Months{t})(U_mask_KDau) = NaN;
    KDau_V_ek.(Months{t})(V_mask_KDau) = NaN;
end


%% UV
KDau_U.mean = KDau_U_g.mean + KDau_U_ek.mean;
KDau_V.mean = KDau_V_g.mean + KDau_V_ek.mean;
KDau_U.mean(U_mask_KDau) = 0;
KDau_V.mean(V_mask_KDau) = 0;
for t = 1 : 12
    KDau_U.(Months{t}) = ...
        KDau_U_g.(Months{t}) + KDau_U_ek.(Months{t});
    KDau_V.(Months{t}) = ...
        KDau_V_g.(Months{t}) + KDau_V_ek.(Months{t});
    KDau_U.(Months{t})(U_mask_KDau) = 0;
    KDau_V.(Months{t})(V_mask_KDau) = 0;
end


%% F
dU.mean = (KDau_U.mean(:,2:end) - KDau_U.mean(:,1:end-1));
for t = 1 : 12
    dU.(Months{t}) = (KDau_U.(Months{t})(:,2:end) - ...
        KDau_U.(Months{t})(:,1:end-1));
end

dx = NaN(size(dU.mean));
for m = 1 : length(lat_u)
    dx_now = a * cos(lat_u(m) * pi180) .* ...
        (lon_u(2:end) - lon_u(1:end-1)) * pi180;
    
    dx(m,:) = dx_now;
end

lat_repmat = repmat(lat_v,1,length(lon_v));

dV.mean = ...
    KDau_V.mean(1:end-1,:) .* ...
    cos(lat_repmat(1:end-1,:) * pi180) - ...
    KDau_V.mean(2:end,:) .* ...
    cos(lat_repmat(2:end,:) * pi180);
for t = 1 : 12
    dV.(Months{t}) = ...
        KDau_V.(Months{t})(1:end-1,:) .* ...
        cos(lat_repmat(1:end-1,:) * pi180) - ...
        KDau_V.(Months{t})(2:end,:) .* ...
        cos(lat_repmat(2:end,:) * pi180);
end

dy = NaN(size(dV.mean));
for jj = 1 : length(lon_v)
    dy_now = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
    
    dy(:,jj) = dy_now;
end
lat_u_repmat = repmat(lat_u,1,length(lon_v));

KDau_F.mean = dU.mean./dx + 1./cos(lat_u_repmat * pi180).*dV.mean./dy;
for t = 1 : 12
    KDau_F.(Months{t}) = dU.(Months{t})./dx + ...
        1./cos(lat_u_repmat * pi180).*dV.(Months{t})./dy;
end

aus8_coor.F_mask_KDau = KDau_F.mean == 0;
save([data_path 'SACS_data/aus8_coor'], 'aus8_coor')


%%
save([data_path 'SACS_data/KDau_U_g_itgr'], 'KDau_U_g')
save([data_path 'SACS_data/KDau_V_g_itgr'], 'KDau_V_g')
save([data_path 'SACS_data/KDau_U_ek'], 'KDau_U_ek')
save([data_path 'SACS_data/KDau_V_ek'], 'KDau_V_ek')
save([data_path 'SACS_data/KDau_U'], 'KDau_U')
save([data_path 'SACS_data/KDau_V'], 'KDau_V')
save([data_path 'SACS_data/KDau_F'], 'KDau_F')

disp(['KDau_U_g_itgr KDau_V_g_itgr KDau_U_ek KDau_V_ek ' ...
    'KDau_U KDau_V KDau_F DONE'])

