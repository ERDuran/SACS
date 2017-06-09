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
for m = 1 : length(lat_u)
    for n = 1 : length(lon_u)
        u_g_now = squeeze(aus8_u_g.mean(m,n,:));
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


% V_g
% integrate V along z. V_z is in m^2/s
for m = 1 : length(lat_v)
    for n = 1 : length(lon_v)
        v_g_now = squeeze(aus8_v_g.mean(m,n,:));
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
% 1) figure set-up
close all
font_size = 8;
fig1 = figure(1);
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
rowN = 4; colN = 2;
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = 0.04; % gap width between subplots
gap_h = 0.04; % gap height between subplots
marg_b = 0.03; % bottom margin
marg_t = 0.03; % top margin
marg_l = 0.02; % left margin
marg_r = 0.02; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length
h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

% 2) all data set-up
sp = 1;
axes(h_axes_sp(sp))
% SBC_U_prime_now = SBC_U_prime;
% SBC_U_prime_now(SBC_U_prime_now<0) = 0;

% 3) asal pcolor set-up
% new colormap routine !
% levels = 11;
% Reds = othercolor('Reds9', levels);
% Reds(1,:) = [1 1 1];
% magnif = 100;
% Reds_cont = ...
%     [0 0.05 0.1 0.5 1 5 10 20 50 100 200 300]*magnif;
% Reds_cont_length = length(Reds_cont);
% cmap_custom = cmapcust(Reds,Reds_cont);
colormap(h_axes_sp(sp), othercolor('RdBu11', 20));

% 4) plot asal pcolor
pcolor(...
    lon_u, ...
    lat_u, ...
    aus8_U_g.mean)
shading interp
% hold on
caxis([-300 300]);
% axis([114 148 -45 -33])
% freezeColors

% cmap = colormap(Reds);
% Reds_linspace = ...
%     linspace(Reds_cont(1), Reds_cont(end), Reds_cont_length);
cbar = colorbar;
% set(cbar, 'YTick',Reds_linspace, 'YTickLabel',Reds_cont/magnif);

% 5) plot SBC boundaries
% hold on
% plot([SC_lon_u_repelem(1,1) SC_lon_u_repelem(1,1)], ...
%     [SC_lat_v_north(1,1) SC_lat_v_south(1,1)], ...
%     'g','linewidth',0.5,'linestyle','-')
% plot(SC_lon_u_repelem(1,:), ...
%     SC_lat_v_north_repelem(1,:), ...
%     'g','linewidth',0.5,'linestyle','-')
% plot(SC_lon_u_repelem(1,:), ...
%     SC_lat_v_south_repelem(1,:), ...
%     'g','linewidth',0.5,'linestyle','-')
% plot([SC_lon_u_repelem(1,end) SC_lon_u_repelem(1,end)], ...
%     [SC_lat_v_north(1,end) SC_lat_v_south(1,end)], ...
%     'g','linewidth',0.5,'linestyle','-')

% 7) UV quiver set-up
% [lon_mg, lat_mg]=meshgrid(lon,lat);
% quiv_S = 5;
% nn = 1;

% 8) plot directional U and V
% quiver(...
%     lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
%     SBC_U_prime_interp2(1:nn:end, 1:nn:end), ...
%     SBC_V_prime_interp2(1:nn:end, 1:nn:end), ...
%     quiv_S, 'k');
% reference arrow goes here

% 9) title, grid, background and fonts
% title('SBC current')
% grid
% set(gca,'xtick',115:2:147)
% set(gca,'layer','top','color',[0.7 0.7 0.7],...
%     'fontsize',font_size,'fontweight','bold')
% if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
% if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end


% % Save
% outputls = ls(figures_path);
% scriptname = mfilename;
% if ~contains(outputls, scriptname)
%     mkdir(figures_path, scriptname)
% end
% export_fig(fig1, [figures_path mfilename '/' scriptname(1:3) '_'], ...
%     '-m4')
% close



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

