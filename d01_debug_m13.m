%%
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_currents'])
load([data_path 'SACS_data/aus8_U_prime'])
load([data_path 'SACS_data/aus8_V_prime'])

a = aus8_coor.a;
pi180 = aus8_coor.pi180;
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;
U_prime_up = aus8_currents.ptop_to_pmid.U_prime.mean;
V_prime_up = aus8_currents.ptop_to_pmid.V_prime.mean;
U_g_prime_dw = aus8_currents.pmid_to_pbot.U_g_prime.mean;
V_g_prime_dw = aus8_currents.pmid_to_pbot.V_g_prime.mean;
U_prime = U_prime_up + U_g_prime_dw;
V_prime = V_prime_up + V_g_prime_dw;
lat_v_SBC_north = aus8_currents.lat_v_SBC_north;
lat_v_SBC_south = aus8_currents.lat_v_SBC_south;
lat_v_DRC_north = aus8_currents.lat_v_DRC_north;
lat_v_DRC_south = aus8_currents.lat_v_DRC_south;
lon_u_ALLC = aus8_currents.lon_u_ALLC;


%% Prep
dx_v = NaN(length(lat_v), length(lon_u)-1);
for ii = 1 : length(lat_v)
    dx_v(ii,:) = a * cos(lat_v(ii) * pi180) * ...
        (lon_u(2:end) - lon_u(1:end-1)) * pi180;
end
dy_raw = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
dy_u = repmat(dy_raw, [1 length(lon_u)]);

%
U_prime_up(isnan(U_prime_up)) = 0;
V_prime_up(isnan(V_prime_up)) = 0;
U_g_prime_dw(isnan(U_g_prime_dw)) = 0;
V_g_prime_dw(isnan(V_g_prime_dw)) = 0;
Ut_prime_up = U_prime_up .* dy_u;
Vt_prime_up = V_prime_up .* dx_v;
Ut_g_prime_dw = U_g_prime_dw .* dy_u;
Vt_g_prime_dw = V_g_prime_dw .* dx_v;
Ut_prime = U_prime .* dy_u;
Vt_prime = V_prime .* dx_v;
%
dx_u = NaN(length(lat_u), length(lon_u)-1);
for ii = 1 : length(lat_u)
    dx_u(ii,:) = a * cos(lat_u(ii) * pi180) * ...
        (lon_u(2:end) - lon_u(1:end-1)) * pi180;
end
dy_raw = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
dy_v = repmat(dy_raw, [1 length(lon_v)]);


%% THIS IS OK !!
aus8_U_prime.mean(isnan(aus8_U_prime.mean)) = 0;
aus8_V_prime.mean(isnan(aus8_V_prime.mean)) = 0;
du = aus8_U_prime.mean(:,2:end) - aus8_U_prime.mean(:,1:end-1);
dudx = du ./ dx_u;
lat_v_repmat = repmat(lat_v,1,length(lon_v));
dv = aus8_V_prime.mean(1:end-1,:).* ...
    cos(lat_v_repmat(1:end-1,:) * pi180) - ...
    aus8_V_prime.mean(2:end,:).* ...
    cos(lat_v_repmat(2:end,:) * pi180);
lat_u_repmat = repmat(lat_u,1,length(lon_v));
dvdy = 1./cos(lat_u_repmat * pi180).*dv./dy_v;
div_UV_prime.mean = (dudx + dvdy) .* dy_v .* dx_u;
div_UV_prime.mean(div_UV_prime.mean==0)=NaN;


du = U_prime_up(:,2:end) - U_prime_up(:,1:end-1);
dudx = du ./ dx_u;
lat_v_repmat = repmat(lat_v,1,length(lon_v));
dv = V_prime_up(1:end-1,:).* ...
    cos(lat_v_repmat(1:end-1,:) * pi180) - ...
    V_prime_up(2:end,:).* ...
    cos(lat_v_repmat(2:end,:) * pi180);
lat_u_repmat = repmat(lat_u,1,length(lon_v));
dvdy = 1./cos(lat_u_repmat * pi180).*dv./dy_v;
div_UV_prime2.mean = (dudx + dvdy) .* dy_v .* dx_u;
div_UV_prime2.mean(div_UV_prime2.mean==0)=NaN;


du = U_g_prime_dw(:,2:end) - U_g_prime_dw(:,1:end-1);
dudx = du ./ dx_u;
lat_v_repmat = repmat(lat_v,1,length(lon_v));
dv = V_g_prime_dw(1:end-1,:).* ...
    cos(lat_v_repmat(1:end-1,:) * pi180) - ...
    V_g_prime_dw(2:end,:).* ...
    cos(lat_v_repmat(2:end,:) * pi180);
lat_u_repmat = repmat(lat_u,1,length(lon_v));
dvdy = 1./cos(lat_u_repmat * pi180).*dv./dy_v;
div_UV_prime3.mean = (dudx + dvdy) .* dy_v .* dx_u;
div_UV_prime3.mean(div_UV_prime3.mean==0)=NaN;

div_UV_prime4.mean = div_UV_prime2.mean + div_UV_prime3.mean;


%%

DIV1 = ...
    Ut_prime(:,2:end) - Ut_prime(:,1:end-1) + ...
    Vt_prime(1:end-1,:) - Vt_prime(2:end,:);

DIV2 = ...
    Ut_prime_up(:,2:end) - Ut_prime_up(:,1:end-1) + ...
    Vt_prime_up(1:end-1,:) - Vt_prime_up(2:end,:);

DIV3 = ...
    Ut_g_prime_dw(:,2:end) - Ut_g_prime_dw(:,1:end-1) + ...
    Vt_g_prime_dw(1:end-1,:) - Vt_g_prime_dw(2:end,:);

DIV4 = ...
    DIV2 + ...
    DIV3;


%%
% dUprime dVprime
aus8_U_prime.mean(isnan(aus8_U_prime.mean)) = 0;
aus8_V_prime.mean(isnan(aus8_V_prime.mean)) = 0;
du = aus8_U_prime.mean(:,2:end) - aus8_U_prime.mean(:,1:end-1);
dudx = du ./ dx_u;
lat_v_repmat = repmat(lat_v,1,length(lon_v));
dv = aus8_V_prime.mean(1:end-1,:).* ...
    cos(lat_v_repmat(1:end-1,:) * pi180) - ...
    aus8_V_prime.mean(2:end,:).* ...
    cos(lat_v_repmat(2:end,:) * pi180);
lat_u_repmat = repmat(lat_u,1,length(lon_v));
dvdy = 1./cos(lat_u_repmat * pi180).*dv./dy_v;
div_UV_prime.mean = dudx + dvdy;

% dUprime dVprime
aus8_U_prime.mean(isnan(aus8_U_prime.mean)) = 0;
aus8_V_prime.mean(isnan(aus8_V_prime.mean)) = 0;
du = aus8_U_prime.mean(:,2:end) - aus8_U_prime.mean(:,1:end-1);
dudx = du ./ dx_u;
lat_v_repmat = repmat(lat_v,1,length(lon_v));
dv = aus8_V_prime.mean(1:end-1,:).* ...
    cos(lat_v_repmat(1:end-1,:) * pi180) - ...
    aus8_V_prime.mean(2:end,:).* ...
    cos(lat_v_repmat(2:end,:) * pi180);
lat_u_repmat = repmat(lat_u,1,length(lon_v));
dvdy = 1./cos(lat_u_repmat * pi180).*dv./dy_v;
div_UV_prime.mean = dudx + dvdy;

%%

