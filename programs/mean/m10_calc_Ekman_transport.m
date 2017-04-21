%% calculate the Ekman transport
clearvars('-except', '*_path')

load ERA_Interim
load aus8_TEOS10
load aus8_ZD_method
U_g = aus8_ZD_method.U_g;
V_g = aus8_ZD_method.V_g;
lat = ERA_Interim.lat;
lon = ERA_Interim.lon;
lat_u = aus8_ZD_method.lat_u;
lon_u = aus8_ZD_method.lon_u;
lat_v = aus8_ZD_method.lat_v;
lon_v = aus8_ZD_method.lon_v;
iews = ERA_Interim.iews.annual;
inss = ERA_Interim.inss.annual;
pres = aus8_TEOS10.pres;
rho = aus8_TEOS10.rho.mean;


%% 
Tau_x = interp2(lon, lat, iews, lon_v, lat_v);
Tau_y = interp2(lon, lat, inss, lon_u, lat_u);

%
f_x = gsw_f(lat_u);
f_x_repmat = repmat(f_x, 1, length(lon_u));
f_y = gsw_f(lat_v);
f_y_repmat = repmat(f_y, 1, length(lon_v));

%
ekman_depth = 30;
ekman_layer_ind = pres <= ekman_depth;

%
rho_0 = nanmean(nanmean(nanmean(rho(:,:,ekman_layer_ind))));

% Southern Hemisphere !
U_ek = Tau_y./f_x_repmat/rho_0;
V_ek = -Tau_x./f_y_repmat/rho_0;

%
pres_mid_length = length(pres) -1;

%
U_ek(isnan(U_g)) = NaN;
V_ek(isnan(V_g)) = NaN;
U_ekman_repmat = repmat(U_ek,1,1,pres_mid_length);
V_ekman_repmat = repmat(V_ek,1,1,pres_mid_length);
U = U_g + U_ek;
V = V_g + V_ek;

%
aus8_ZD_method.Tau_x = Tau_x;
aus8_ZD_method.Tau_y = Tau_y;
aus8_ZD_method.ekman_depth = ekman_depth;
aus8_ZD_method.U_ek = U_ek;
aus8_ZD_method.V_ek = V_ek;
aus8_ZD_method.U = U;
aus8_ZD_method.V = V;


%%
save([cars_out_path 'aus8_ZD_method'], 'aus8_ZD_method')
disp(['aus8_ZD_method saved in ' ...
    cars_out_path 'aus8_ZD_method'])

