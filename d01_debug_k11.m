%%
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/KDau_currents'])
load([data_path 'SACS_data/aus8_currents'])
load([data_path 'SACS_data/KDau_U_prime'])
load([data_path 'SACS_data/KDau_V_prime'])
load([data_path 'SACS_data/KDau_u_g'])
load([data_path 'SACS_data/KDau_v_g'])
load([data_path 'SACS_data/KDau_u_g_prime'])
load([data_path 'SACS_data/KDau_v_g_prime'])

a = aus8_coor.a;
pi180 = aus8_coor.pi180;
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;
U_prime_up = KDau_currents.ptop_to_pmid.U_prime.mean;
V_prime_up = KDau_currents.ptop_to_pmid.V_prime.mean;
U_g_prime_dw = KDau_currents.pmid_to_pbot.U_g_prime.mean;
V_g_prime_dw = KDau_currents.pmid_to_pbot.V_g_prime.mean;
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


%% CALCULATIONS ARE OK !!
KDau_U_prime.mean(isnan(KDau_U_prime.mean)) = 0;
KDau_V_prime.mean(isnan(KDau_V_prime.mean)) = 0;
du = KDau_U_prime.mean(:,2:end) - KDau_U_prime.mean(:,1:end-1);
dudx = du ./ dx_u;
lat_v_repmat = repmat(lat_v,1,length(lon_v));
dv = KDau_V_prime.mean(1:end-1,:).* ...
    cos(lat_v_repmat(1:end-1,:) * pi180) - ...
    KDau_V_prime.mean(2:end,:).* ...
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

DIV1 = ...
    (Ut_prime(:,2:end) - Ut_prime(:,1:end-1) + ...
    Vt_prime(1:end-1,:) - Vt_prime(2:end,:)) * 10^-6;

DIV2 = ...
    Ut_prime_up(:,2:end) - Ut_prime_up(:,1:end-1) + ...
    Vt_prime_up(1:end-1,:) - Vt_prime_up(2:end,:);

DIV3 = ...
    Ut_g_prime_dw(:,2:end) - Ut_g_prime_dw(:,1:end-1) + ...
    Vt_g_prime_dw(1:end-1,:) - Vt_g_prime_dw(2:end,:);

DIV4 = ...
    DIV2 + ...
    DIV3;


%% START OF FC... ERROR NOT IN U
lat_v_DRC_north_1 = find(lat_v == lat_v_DRC_north(end));

lon_u_DRC_1 = find(lon_u == 147);
lon_u_DRC_2 = find(lon_u == 146.875);

lat_v_DRC_south_1 = find(lat_v == lat_v_DRC_south(end))-1;

lat_v_SBC_north_1 = find(lat_v == lat_v_SBC_north(end));

FC_U_up1 = Ut_prime_up(lat_v_DRC_north_1:lat_v_DRC_south_1,lon_u_DRC_1);
FC_U_up1_trans = sum(FC_U_up1) * 10^-6;

FC_U_dw1 = Ut_g_prime_dw(lat_v_SBC_north_1:lat_v_DRC_south_1,lon_u_DRC_1);
FC_U_dw1_trans = sum(FC_U_dw1) * 10^-6;

FC_U_up2 = Ut_prime_up(lat_v_DRC_north_1:lat_v_DRC_south_1,lon_u_DRC_2);
FC_U_up2_trans = sum(FC_U_up2) * 10^-6;

FC_U_dw2 = Ut_g_prime_dw(lat_v_SBC_north_1:lat_v_DRC_south_1,lon_u_DRC_2);
FC_U_dw2_trans = sum(FC_U_dw2) * 10^-6;

FC_U_1 = FC_U_up1_trans + FC_U_dw1_trans;
FC_U_2 = FC_U_up2_trans + FC_U_dw2_trans;


%% START OF FC... ERROR NOT IN UPPER V

FC_Vn_up1 = -Vt_prime_up(lat_v_DRC_north_1,lon_u_DRC_2);
FC_Vn_up1_trans = sum(FC_Vn_up1) * 10^-6;

FC_Vs_up1 = Vt_prime_up(lat_v_DRC_south_1+1,lon_u_DRC_2);
FC_Vs_up1_trans = sum(FC_Vs_up1) * 10^-6;


%% START OF FC... ERROR NOT IN LOWER V

FC_Vn_dw1 = -Vt_g_prime_dw(lat_v_SBC_north_1,lon_u_DRC_2);
FC_Vn_dw1_trans = sum(FC_Vn_dw1) * 10^-6;

FC_Vs_dw1 = Vt_g_prime_dw(lat_v_DRC_south_1+1,lon_u_DRC_2);
FC_Vs_dw1_trans = sum(FC_Vs_dw1) * 10^-6;


%% START OF FC... ERROR NOT IN V !
FC_Vn = FC_Vn_up1_trans + FC_Vn_dw1_trans;
FC_Vs = FC_Vs_up1_trans + FC_Vs_dw1_trans;


%% START OF FC... ERROR NOT IN W
nc = 0;
nc = nc + 1;
lat_v_SBC_ind = lat_v_SBC_north_1 : lat_v_DRC_north_1;

lat_u_SBC_north_ind = find(lat_u<lat_v_SBC_north(end), 1, 'first');
lat_u_SBC_south_ind = find(lat_u>lat_v_DRC_north(end), 1, 'last');
lat_u_SBC_ind = lat_u_SBC_north_ind : lat_u_SBC_south_ind;

V_prime_up_SBC_for_W = ...
    V_prime_up(lat_v_SBC_ind,lon_u_DRC_2);
U_prime_up_SBC_for_W = ...
    U_prime_up(lat_u_SBC_ind,lon_u_DRC_2:lon_u_DRC_1);

du = U_prime_up_SBC_for_W(:,2:end) - ...
    U_prime_up_SBC_for_W(:,1:end-1);
dudx = du ./ dx_u(lat_u_SBC_ind,lon_u_DRC_2);
dv = V_prime_up_SBC_for_W(1:end-1,:).* ...
    cos(lat_v_repmat(lat_u_SBC_ind,lon_u_DRC_2) * pi180) - ...
    V_prime_up_SBC_for_W(2:end,:).* ...
    cos(lat_v_repmat(lat_u_SBC_ind+1,lon_u_DRC_2) * pi180);
dvdy = 1./cos(lat_u_repmat(lat_u_SBC_ind,lon_u_DRC_2) * pi180).*dv./...
    dy_v(lat_u_SBC_ind,lon_u_DRC_2);
div_UV_prime_up_SBC = dudx + dvdy;

% for ii = 1 : length(lat_u)
%     for jj = 1 : length(lon_v)
%         if isfinite(div_UV_prime_up_SBC(ii,jj))
%             if div_UV_prime_up_SBC(ii,jj) == div_UV_prime.mean(ii,jj)
%                 div_UV_prime_up_SBC(ii,jj) = 0;
%             end
%         end
%     end
% end
% disp('div_UV_prime_up_SBC was corrected')

w_prime_up_SBC = div_UV_prime_up_SBC;
W_prime_up_mcps_SBC = w_prime_up_SBC .* ...
    dy_v(lat_u_SBC_ind,lon_u_DRC_2) .* dx_u(lat_u_SBC_ind,lon_u_DRC_2);
SBC_Wt_mcps = nansum(W_prime_up_mcps_SBC,1);
SBC_Wt = SBC_Wt_mcps * 10^-6;


%%
nc = 0;
nc = nc + 1;
lat_v_SBC_ind = lat_v_SBC_north_1 : lat_v_DRC_north_1;

lat_u_SBC_north_ind = find(lat_u<lat_v_SBC_north(end), 1, 'first');
lat_u_SBC_south_ind = find(lat_u>lat_v_DRC_north(end), 1, 'last');
lat_u_SBC_ind = lat_u_SBC_north_ind : lat_u_SBC_south_ind;

Vt_prime_up_SBC_for_W = ...
    Vt_prime_up(lat_v_SBC_ind,lon_u_DRC_2);
Ut_prime_up_SBC_for_W = ...
    Ut_prime_up(lat_u_SBC_ind,lon_u_DRC_2:lon_u_DRC_1);

du = Ut_prime_up_SBC_for_W(:,2:end) - ...
    Ut_prime_up_SBC_for_W(:,1:end-1);
dv = Vt_prime_up_SBC_for_W(1:end-1,:) - ...
    Vt_prime_up_SBC_for_W(2:end,:);

w_prime_up_SBC = du + dv;
% W_prime_up_mcps_SBC = w_prime_up_SBC .* ...
%     dy_v(lat_u_SBC_ind,lon_u_DRC_2) .* dx_u(lat_u_SBC_ind,lon_u_DRC_2);
SBC_Wt_mcps = nansum(w_prime_up_SBC,1);
SBC_Wt2 = SBC_Wt_mcps * 10^-6;


%% START OF FC... ERROR NOT IN W
% FC

FC_star = FC_U_1 -(-SBC_Wt + FC_Vn + FC_Vs)


%%
nc = 0;
nc = nc + 1;
lat_v_SBC_ind = lat_v_SBC_north_1 : lat_v_DRC_north_1;

lat_u_SBC_north_ind = find(lat_u<lat_v_SBC_north(end), 1, 'first');
lat_u_SBC_south_ind = find(lat_u>lat_v_DRC_north(end), 1, 'last');
lat_u_SBC_ind = lat_u_SBC_north_ind : lat_u_SBC_south_ind;

Vt_g_prime_dw_SBC_for_W = ...
    Vt_g_prime_dw(lat_v_SBC_ind,lon_u_DRC_2);
Ut_g_prime_dw_SBC_for_W = ...
    Ut_g_prime_dw(lat_u_SBC_ind,lon_u_DRC_2:lon_u_DRC_1);

du = Ut_g_prime_dw_SBC_for_W(:,2:end) - ...
    Ut_g_prime_dw_SBC_for_W(:,1:end-1);
dv = Vt_g_prime_dw_SBC_for_W(1:end-1,:) - ...
    Vt_g_prime_dw_SBC_for_W(2:end,:);

w_g_prime_dw_SBC = du + dv;
% W_prime_up_mcps_SBC = w_prime_up_SBC .* ...
%     dy_v(lat_u_SBC_ind,lon_u_DRC_2) .* dx_u(lat_u_SBC_ind,lon_u_DRC_2);
SBC_Wt_mcps = nansum(w_g_prime_dw_SBC,1);
SBC_Wt3 = SBC_Wt_mcps * 10^-6;


FC_star = FC_U_1 -(SBC_Wt3 + FC_Vn + FC_Vs)

