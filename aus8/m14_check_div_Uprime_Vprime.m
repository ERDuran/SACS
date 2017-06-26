%% check that div(U_zd,V_zd) = 0
clearvars('-except', '*_path')

load aus8_ZD_method
depth_thicknesses = aus8_ZD_method.depth_thicknesses;
u_g_prime = aus8_ZD_method.u_g_prime;
v_g_prime = aus8_ZD_method.v_g_prime;
a = aus8_ZD_method.a;
pi180 = aus8_ZD_method.pi180;
n_it = aus8_ZD_method.n_iteration;
lat_u = aus8_ZD_method.lat_u;
lon_u = aus8_ZD_method.lon_u;
lat_v = aus8_ZD_method.lat_v;
lon_v = aus8_ZD_method.lon_v;
U_ek = aus8_ZD_method.U_ek;
V_ek = aus8_ZD_method.V_ek;
F_nan = aus8_ZD_method.F == 0;


%% calculate vertical integration
lat_u_length = length(lat_u);
lon_u_length = length(lon_u);
U_g_prime = NaN(lat_u_length, lon_u_length);
% integrate U along z. U_z is in m^2/s
for ii = 1 : lat_u_length
    for jj = 1 : lon_u_length
        u_mid_now = squeeze(u_g_prime(ii,jj,:));
        u_mid_now_times_depth_thicknesses = ...
            u_mid_now .* depth_thicknesses;
        U_g_prime(ii,jj) = nansum(u_mid_now_times_depth_thicknesses);
    end
end
u_mid_mask = isnan(u_g_prime(:,:,1));
U_g_prime(u_mid_mask) = NaN;

lat_v_length = length(lat_v);
lon_v_length = length(lon_v);
V_g_prime = NaN(lat_v_length, lon_v_length);
% integrate V along z. V_z is in m^2/s
for ii = 1 : lat_v_length
    for jj = 1 : lon_v_length
        v_mid_now = squeeze(v_g_prime(ii,jj,:));
        v_mid_now_times_depth_thicknesses = ...
            v_mid_now .* depth_thicknesses;
        V_g_prime(ii,jj) = nansum(v_mid_now_times_depth_thicknesses);
    end
end
v_mid_mask = isnan(v_g_prime(:,:,1));
V_g_prime(v_mid_mask) = NaN;

% Re-add the Ekman drift
U_prime = U_g_prime + U_ek;
V_prime = V_g_prime + V_ek;
U_prime(isnan(U_g_prime)) = 0;
V_prime(isnan(V_g_prime)) = 0;

% calculate divergence
du = (U_prime(:,2:end) - U_prime(:,1:end-1));
dx = NaN(size(du));
for ii = 1 : length(lat_u)
    dx_now = a * cos(lat_u(ii) * pi180) .* ...
        (lon_u(2:end) - lon_u(1:end-1)) * pi180;
    dx(ii,:) = dx_now;
end
lat_v_repmat = repmat(lat_v,1,length(lon_v));
dv = V_prime(1:end-1,:).*cos(lat_v_repmat(1:end-1,:) * pi180) - ...
    V_prime(2:end,:).*cos(lat_v_repmat(2:end,:) * pi180);
dy = NaN(size(dv));
for jj = 1 : length(lon_v)
    dy_now = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
    dy(:,jj) = dy_now;
end
lat_u_repmat = repmat(lat_u,1,length(lon_v));
div_UV_prime = du./dx + 1./cos(lat_u_repmat * pi180).*dv./dy;
U_prime(isnan(U_g_prime)) = NaN;
V_prime(isnan(V_g_prime)) = NaN;
div_UV_prime(F_nan) = NaN;


%% Plot setup
lat_div = aus8_ZD_method.lat_F;
lon_div = aus8_ZD_method.lon_F;
F = aus8_ZD_method.F;
F(F == 0) = NaN;
Lap_phi = aus8_ZD_method.Lap_phi;
lat_u = aus8_ZD_method.lat_u;
lon_u = aus8_ZD_method.lon_u;
U = aus8_ZD_method.U;
U_d = aus8_ZD_method.U_d;
U_d(U_d == 0) = NaN;
lat_v = aus8_ZD_method.lat_v;
lon_v = aus8_ZD_method.lon_v;
V = aus8_ZD_method.V;
V_d = aus8_ZD_method.V_d;
V_d(V_d == 0) = NaN;
diff_F__Lap_phi = F - Lap_phi;
U_interp2 = interp2(lon_u,lat_u,U,lon_div,lat_div);
V_interp2 = interp2(lon_v,lat_v,V,lon_div,lat_div);
U_d_interp2 = interp2(lon_u,lat_u,U_d,lon_div,lat_div);
V_d_interp2 = interp2(lon_v,lat_v,V_d,lon_div,lat_div);
U_prime_interp2 = interp2(lon_u,lat_u,U_prime,lon_div,lat_div);
V_prime_interp2 = interp2(lon_v,lat_v,V_prime,lon_div,lat_div);
UV_speed = sqrt(U_interp2.^2 + V_interp2.^2);
UV_d_speed = sqrt(U_d_interp2.^2 + V_d_interp2.^2);
UV_prime_speed = sqrt(U_prime_interp2.^2 + V_prime_interp2.^2);


%% plot map one order
close
fig1 = figure;
x_pos = 0; % x figure position
y_pos = 0.05; % y figure position
x_lgth = 0.95; % x figure length
y_lgth = 0.99; % y figure length
set(gcf,'units','normalized','outerposition', ...
    [x_pos y_pos x_lgth y_lgth]) % apply onto current figure

rowN = 2; % number of rows
colN = 3; % number of columns
gap_w = 0.02; % gap width between subplots
gap_h = 0.05; % gap height between subplots
marg_b = 0.05; % bottom margin
marg_t = 0.04; % top margin
marg_l = 0.025; % left margin
marg_r = 0.01; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length
h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]); % creates subplots

font_size = 8; % font size
background_c = [0.7 0.7 0.7]; % background colour

%
ax_now = h_axes_sp(1); axes(ax_now) % call one subplot;
order = 5; % cbar order of magnitude
pcolor(lon_div,lat_div,F), shading interp
caxis([-10^-order 10^-order])
n_levs = 20;
colormap(ax_now, othercolor('RdGy11', n_levs))
colorbar
set(gca,'color',background_c)
title(sprintf(['F']))
grid
set(gca,'layer','top', 'tickdir', 'out', ...
    'fontsize',font_size, 'fontweight','bold')
set(gca,'xticklabel','')
% set(gca,'yticklabel','')

%
ax_now = h_axes_sp(2); axes(ax_now) % call one subplot;
order = 6; % cbar order of magnitude
pcolor(lon_div,lat_div,diff_F__Lap_phi), shading interp
caxis([-10^-order 10^-order])
n_levs = 20;
colormap(ax_now, othercolor('RdGy11', n_levs))
colorbar
set(gca,'color',background_c)
title(sprintf(['F - Lap(phi)']))
grid
set(gca,'layer','top', 'tickdir', 'out', ...
    'fontsize',font_size, 'fontweight','bold')
set(gca,'xticklabel','')
set(gca,'yticklabel','')

%
ax_now = h_axes_sp(3); axes(ax_now) % call one subplot;
pcolor(lon_div,lat_div,div_UV_prime), shading interp
caxis([-10^-order 10^-order])
n_levs = 20;
colormap(ax_now, othercolor('RdGy11', n_levs))
colorbar
set(gca,'color',background_c)
title(sprintf(['div(Uprime,Vprime)']))
grid
set(gca,'layer','top', 'tickdir', 'out', ...
    'fontsize',font_size, 'fontweight','bold')
set(gca,'xticklabel','')
set(gca,'yticklabel','')

%
ax_now = h_axes_sp(4); axes(ax_now) % call one subplot;
pcolor(lon_div,lat_div,UV_speed), shading interp
order = 30;
caxis([0 order])
n_levs = 12;
colormap(ax_now, flipud(othercolor('Blues9', n_levs)))
colorbar
hold on
nn = 4;
quiv_S = 4;
[lon_mg, lat_mg]=meshgrid(lon_div,lat_div);
h_quiv = quiver(...
    lon_mg(1:nn:end,1:nn:end), lat_mg(1:nn:end,1:nn:end), ...
    U_interp2(1:nn:end,1:nn:end), ...
    V_interp2(1:nn:end,1:nn:end), ...
    quiv_S, 'k');
set(gca,'color',background_c)
title(sprintf(['U,V']))
grid
set(gca,'layer','top', 'tickdir', 'out', ...
    'fontsize',font_size, 'fontweight','bold')
% set(gca,'xticklabel','')
% set(gca,'yticklabel','')
gca_pos = get(gca,'position');
h = axes('Position', ...
    gca_pos+[x_sp-0.114 y_sp-0.082 -x_sp+0.096 -y_sp+0.041],...
    'Layer','top');
[u_ref, v_ref] = deal(NaN(5,5));
u_ref(3,3) = max(max(UV_speed))/2;
v_ref(3,3) = 0;
u_ref(15,15) = max(max(UV_speed));
v_ref(15,15) = 0;
hquiv = quiver(...
    lon_mg(1:nn:nn*14+1,1:nn:nn*14+1)+1.5, ...
    lat_mg(1:nn:nn*14+1,1:nn:nn*14+1)+1, ...
    u_ref, ...
    v_ref, ...
    quiv_S,'k');
text(lon_div(10),lat_div(4)-0.06,[num2str(round(u_ref(3,3)*100)/100) ...
    ' m^2/s'], ...
    'fontsize',font_size)
axis([108 118 -31.5 -29.5])
set(h,'xticklabel','','yticklabel','')

%
ax_now = h_axes_sp(5); axes(ax_now) % call one subplot;
pcolor(lon_div,lat_div,UV_d_speed), shading interp
caxis([0 order])
colormap(ax_now, flipud(othercolor('Blues9', n_levs)))
colorbar
hold on
[lon_mg, lat_mg]=meshgrid(lon_div,lat_div);
h_quiv = quiver(...
    lon_mg(1:nn:end,1:nn:end), lat_mg(1:nn:end,1:nn:end), ...
    U_d_interp2(1:nn:end,1:nn:end), ...
    V_d_interp2(1:nn:end,1:nn:end), ...
    quiv_S, 'k');
set(gca,'color',background_c)
title(sprintf(['Ud,Vd']))
grid
set(gca,'layer','top', 'tickdir', 'out', ...
    'fontsize',font_size, 'fontweight','bold')
% set(gca,'xticklabel','')
set(gca,'yticklabel','')
gca_pos = get(gca,'position');
h = axes('Position', ...
    gca_pos+[x_sp-0.114 y_sp-0.082 -x_sp+0.096 -y_sp+0.041],...
    'Layer','top');
[u_ref, v_ref] = deal(NaN(5,5));
u_ref(3,3) = max(max(UV_d_speed))/2;
v_ref(3,3) = 0;
u_ref(15,15) = max(max(UV_d_speed));
v_ref(15,15) = 0;
hquiv = quiver(...
    lon_mg(1:nn:nn*14+1,1:nn:nn*14+1)+1.5, ...
    lat_mg(1:nn:nn*14+1,1:nn:nn*14+1)+1, ...
    u_ref, ...
    v_ref, ...
    quiv_S,'k');
text(lon_div(10),lat_div(4)-0.06,[num2str(round(u_ref(3,3)*100)/100) ...
    ' m^2/s'], ...
    'fontsize',font_size)
axis([108 118 -31.5 -29.5])
set(h,'xticklabel','','yticklabel','')

%
ax_now = h_axes_sp(6); axes(ax_now) % call one subplot;
pcolor(lon_div,lat_div,UV_prime_speed), shading interp
caxis([0 order])
colormap(ax_now, flipud(othercolor('Blues9', n_levs)))
colorbar
hold on
[lon_mg, lat_mg]=meshgrid(lon_div,lat_div);
h_quiv = quiver(...
    lon_mg(1:nn:end,1:nn:end), lat_mg(1:nn:end,1:nn:end), ...
    U_prime_interp2(1:nn:end,1:nn:end), ...
    V_prime_interp2(1:nn:end,1:nn:end), ...
    quiv_S, 'k');
set(gca,'color',background_c)
title(sprintf(['Uprime,Vprime']))
grid
set(gca,'layer','top', 'tickdir', 'out', ...
    'fontsize',font_size, 'fontweight','bold')
% set(gca,'xticklabel','')
set(gca,'yticklabel','')
gca_pos = get(gca,'position');
h = axes('Position', ...
    gca_pos+[x_sp-0.114 y_sp-0.082 -x_sp+0.096 -y_sp+0.041],...
    'Layer','top');
[u_ref, v_ref] = deal(NaN(5,5));
u_ref(3,3) = max(max(UV_prime_speed))/2;
v_ref(3,3) = 0;
u_ref(15,15) = max(max(UV_prime_speed));
v_ref(15,15) = 0;
hquiv = quiver(...
    lon_mg(1:nn:nn*14+1,1:nn:nn*14+1)+1.5, ...
    lat_mg(1:nn:nn*14+1,1:nn:nn*14+1)+1, ...
    u_ref, ...
    v_ref, ...
    quiv_S,'k');
text(lon_div(10),lat_div(4)-0.06,[num2str(round(u_ref(3,3)*100)/100) ...
    ' m^2/s'], ...
    'fontsize',font_size)
axis([108 118 -31.5 -29.5])
set(h,'xticklabel','','yticklabel','')

outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
export_fig(fig1, [figures_path mfilename '/' scriptname(1:3) '_'], ...
    '-m4')
close


%%
aus8_ZD_method.U_prime = U_prime;
aus8_ZD_method.V_prime = V_prime;
aus8_ZD_method.div_UV_prime = div_UV_prime;

save([cars_out_path 'aus8_ZD_method'], 'aus8_ZD_method')
disp(['aus8_ZD_method saved in ' ...
    cars_out_path 'aus8_ZD_method'])

