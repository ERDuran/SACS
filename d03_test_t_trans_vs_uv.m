%%
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])

lat = aus8_coor.lat;
lon = aus8_coor.lon;
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;
depth_mid = aus8_coor.depth_mid;
depth = aus8_coor.depth;
a = aus8_coor.a;
pi180 = aus8_coor.pi180;
MTH = aus8_coor.MTH;


%%
for t = 1 : 4
    [data.(MTH{t}), ~, kds_att] = ...
        nc2mat([data_path ...
        'KDS75/KDS75_ncra_y103to109_' MTH{t} '.nc'], ...
        'ALL');
    
    disp([MTH{t} ' OK!'])
end
% fix axis
yt_ocean = flipud(data.JFM.yt_ocean);
xt_ocean = (data.JFM.xt_ocean + 360)';
st_ocean = -data.JFM.st_ocean;
yu_ocean = flipud(data.JFM.yu_ocean);
xu_ocean = (data.JFM.xu_ocean + 360)';
sw_ocean = -data.JFM.sw_ocean;

% u, st_ocean, yu_ocean, xu_ocean
% v, st_ocean, yu_ocean, xu_ocean
% wt, sw_ocean, yt_ocean, xt_ocean
% tx_trans, st_ocean, yt_ocean, xu_ocean
% ty_trans, st_ocean, yu_ocean, xt_ocean

u = flipud(permute(data.JFM.u, [2, 1, 3]));
v = flipud(permute(data.JFM.v, [2, 1, 3]));
w = flipud(permute(data.JFM.wt, [2, 1, 3]));
tx_trans = flipud(permute(data.JFM.tx_trans, [2, 1, 3]));
ty_trans = flipud(permute(data.JFM.ty_trans, [2, 1, 3]));


%%
m = 20;
n = 20;

u_mid = u(m-1,n+1,1)
tx_trans_mid = tx_trans(m,n+1,1)

dz = -st_ocean(1)*2;
dy = a * (yu_ocean(1) - yu_ocean(2)) * pi180;

u_mid_trans = u_mid * dy * dz * 1035 * 10^-9
tx = tx_trans_mid / dy / dz / 10^-6


%%
n = 1;

% u_mid = data.JFM.u(1,1,1) + (data.JFM.u(1,1,2) - data.JFM.u(1,1,1))
v_mid = v(1,n+1,1)
ty_trans_mid = ty_trans(n,1,1)

dz = -st_ocean(1)*2;
dx = a * cos(yt_ocean(1) * pi180) * (xt_ocean(2) - xt_ocean(1)) * pi180;
    
v_mid_trans = v_mid * dx * dz * 10^-6
ty = ty_trans_mid / dx / dz / 10^-6



%%
m = 20;
n = 20;

u_mid = u(m-1,n+1,1) - u(m-1,n,1)

v_mid = v(m-1,n+1,1) - v(m,n+1,1)

w_mid = w(m,n,1)

dz = -st_ocean(1)*2;
dy = a * (yu_ocean(m) - yu_ocean(m+1)) * pi180;
dx = a * cos(yt_ocean(m) * pi180) * (xt_ocean(n+1) - xt_ocean(n)) * pi180;
    
dudx = u_mid * dx;
dvdy = v_mid * dx;


%%
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

% aus8_coor.F_mask_KDau = KDau_F.mean == 0;
% save([data_path 'SACS_data/aus8_coor'], 'aus8_coor')



