%% define all levels of current boundaries
clearvars('-except', 'outputpath')
load aus8_ZD_method
load aus8_currents


%% load data
% note since I want to plot a map using quiver,
% u and v must be on the same grid (T and S grid)
% therefore this will ruin the zero divergence property
% which is ok here because we just want to look at maps of
% geostrophic velocity fields

pres = (0:10:2000)';
pres_mid = aus8_ZD_method.pres_mid;

lat_u = aus8_ZD_method.lat_u;
lon_u = aus8_ZD_method.lon_u;
u_g_prime = aus8_ZD_method.u_g_prime;

lat_v = aus8_ZD_method.lat_v;
lon_v = aus8_ZD_method.lon_v;
v_g_prime = aus8_ZD_method.v_g_prime;

lat_F = aus8_ZD_method.lat_F;
lon_F = aus8_ZD_method.lon_F;
div_UV_prime = aus8_ZD_method.div_UV_prime;

ALLC_west_lon_u = aus8_currents.ALLC_west_lon_u;
ALLC_east_lon_u = aus8_currents.ALLC_east_lon_u;
ALLC_west_to_east_lon_u_ind = aus8_currents.ALLC_west_to_east_lon_u_ind;

SC_bottom_pres = aus8_currents.SC.bottom.pres;
% SC_bottom_lat_u = aus8_currents.SC.bottom.lat_u;
% SC_bottom_lon_u = aus8_currents.SC.bottom.lon_u;
% SC_bottom_lat_v = aus8_currents.SC.bottom.lat_v;
% SC_bottom_lon_v = aus8_currents.SC.bottom.lon_v;
% SC_bottom_lat_v_south_line = aus8_currents.SC.bottom.lat_v_south_line;


%% get all levels above bottom
%
SC_u = zeros(size(u_g_prime));
SC_v = zeros(size(v_g_prime));
SC_pres_mid_length = length(find(pres_mid < SC_bottom_pres));

%
jj_count = 0;
for jj = ALLC_west_to_east_lon_u_ind
    jj_count = jj_count +1;
    
    %
    SC_bottom_lat_u_ind_now = ...
        find(ismember(lat_u, SC_bottom_lat_u(:,jj_count)));
    %
    SC_u(...
        SC_bottom_lat_u_ind_now, ...
        jj, ...
        SC_pres_mid_length(end)) = u_g_prime(...
        SC_bottom_lat_u_ind_now, ...
        jj, ...
        SC_pres_mid_length(end));
    
    if jj ~= ALLC_west_to_east_lon_u_ind(end)
        %
        SC_bottom_lat_v_ind_now = ...
            find(ismember(lat_v, SC_bottom_lat_v(:,jj_count)));
        %
        SC_v(...
            SC_bottom_lat_v_ind_now, ...
            jj, ...
            SC_pres_mid_length(end)) = v_g_prime(...
            SC_bottom_lat_v_ind_now, ...
            jj, ...
            SC_pres_mid_length(end));
    end
end


%%% GET U AND V
for kk = SC_pres_mid_length-1:-1:1
    % longitude index counter
    jj_count = 0;
    % for all longitudes within the west and east boundaries
    for jj = ALLC_west_to_east_lon_u_ind
        % increase lon counter
        jj_count = jj_count + 1;
        
        % lat u
        SC_lat_u_now = u_g_prime(:,jj,kk);
        
        if jj >= 57 && jj <= 65
            SC_lat_u_now(1:32) = NaN;
        end
        
        
        flag = 0;
        jj_1 = 129;
        if jj >= 57 && jj <= jj_1
            p_now = 15; p_now_ind = find(pres_mid == p_now)+1;
            if pres_mid(kk) <= p_now
                flag = 1;
                mask_ind = SC_v(:,jj,p_now_ind) == 0;
                SC_lat_u_now(mask_ind) = NaN;
            end
            
            jj_2 = 138;
        elseif jj > jj_1 && jj <= jj_2
            p_now = 105; p_now_ind = find(pres_mid == p_now)+1;
            if pres_mid(kk) <= p_now
                flag = 1;
                mask_ind = SC_v(:,jj,p_now_ind) == 0;
                SC_lat_u_now(mask_ind) = NaN;
            end
            
        elseif jj > jj_2
            p_now = 105; p_now_ind = find(pres_mid == p_now)+1;
            if pres_mid(kk) <= p_now
                flag = 1;
                mask_ind = SC_v(:,jj,p_now_ind) == 0;
                SC_lat_u_now(mask_ind) = NaN;
            end
        end
        
        SC_lat_u_first_ocean_now_ind = ...
            find(isfinite(SC_lat_u_now), 1, 'first');
        
        
        if jj ~= ALLC_west_to_east_lon_u_ind(end)
            % south
            SC_lat_u_south_ind_now = ...
                find(ismember(lat_u, ...
                SC_bottom_lat_v_south_line(jj_count)+1/16));
            
            SC_lat_u_north_south_vec_ind_now = ...
                (SC_lat_u_first_ocean_now_ind : ...
                SC_lat_u_south_ind_now)';
            
            %
            SC_lat_v_north_south_vec_ind_now = ...
                SC_lat_u_north_south_vec_ind_now;
            SC_lat_v_north_south_vec_ind_now(end+1) = ...
                SC_lat_u_north_south_vec_ind_now(end)+1; %#ok<SAGROW>
            v_g_prime_now = v_g_prime(...
                SC_lat_v_north_south_vec_ind_now, ...
                jj, ...
                kk);
            SC_v(...
                SC_lat_v_north_south_vec_ind_now, ...
                jj, ...
                kk) = v_g_prime_now;
        end
                
        %
        u_g_prime_now = u_g_prime(...
            SC_lat_u_north_south_vec_ind_now, ...
            [jj, jj+1], ...
            kk);
        SC_u(...
            SC_lat_u_north_south_vec_ind_now, ...
            [jj, jj+1], ...
            kk) = u_g_prime_now;
        
    end
end
%%% GET U AND V

SC_u_ind = SC_u ~= 0;
SC_v_ind = SC_v ~= 0;


% GET NORTH SOUTH BOUNDARIES FOR CONTOUR INTEGRATION
% lon u
SC_lon_u_all = repmat(lon_u, [length(lat_u), 1, length(pres_mid)]);
SC_lon_u_all(~SC_u_ind) = NaN;
SC_lon_u = repmat(...
    lon_u(ALLC_west_to_east_lon_u_ind), ...
    [length(pres_mid), 1]);

% lat v
SC_lat_v_all = repmat(lat_v, [1, length(lon_v), length(pres_mid)]);
SC_lat_v_all(~SC_v_ind) = NaN;
[SC_lat_v_north, SC_lat_v_south] = ...
    deal(NaN(length(pres_mid), length(ALLC_west_to_east_lon_u_ind)-1));
for kk = 1 : length(pres_mid)
    jj_count = 0;
    for jj = ALLC_west_to_east_lon_u_ind
        jj_count = jj_count + 1;
        SC_lat_v_north_ind_now = ...
            find(isfinite(SC_lat_v_all(:,jj,kk)), 1, 'first');
        SC_lat_u_south_ind_now = ...
            find(isfinite(SC_lat_v_all(:,jj,kk)), 1, 'last');
        
        if ~isempty(SC_lat_v_north_ind_now)
            SC_lat_v_north(kk, jj_count) = ...
                SC_lat_v_all(SC_lat_v_north_ind_now, jj, kk);
            SC_lat_v_south(kk, jj_count) = ...
                SC_lat_v_all(SC_lat_u_south_ind_now, jj, kk);
        end
    end
end

%
SC_lon_u_repelem = [...
    SC_lon_u(:,1), ...
    repelem(SC_lon_u(:,2:end-1), 1, 2), ...
    SC_lon_u(:,end)];
SC_lat_v_north_repelem = ...
    repelem(SC_lat_v_north, 1, 2);
SC_lat_v_south_repelem = ...
    repelem(SC_lat_v_south, 1, 2);


%%% SAVE
aus8_currents.SC.u_g_prime_ind = SC_u_ind;
aus8_currents.SC.v_g_prime_ind = SC_v_ind;
aus8_currents.SC.lon_u = SC_lon_u;
aus8_currents.SC.lat_v_north = SC_lat_v_north;
aus8_currents.SC.lat_v_south = SC_lat_v_south;
aus8_currents.SC.lon_u_repelem = SC_lon_u_repelem;
aus8_currents.SC.lat_v_north_repelem = SC_lat_v_north_repelem;
aus8_currents.SC.lat_v_south_repelem = SC_lat_v_south_repelem;
%%% SAVE


% check
%
pres_mid_now = 5;
if pres_mid_now == 5
    quiv_S = 10;
else
    quiv_S = 3;
end

p_ind = find(pres_mid == pres_mid_now);

%
lon_plot_west = 114;
lon_plot_east = 149;
lon_F_ind = find(lon_F >= lon_plot_west & lon_F <= lon_plot_east);
lon_u_ind = find(lon_u >= lon_plot_west & lon_u <= lon_plot_east);
lon_v_ind = find(lon_v >= lon_plot_west & lon_v <= lon_plot_east);
lat_plot_north = -31.5;
lat_plot_south = -45;
lat_F_ind = find(lat_F <= lat_plot_north & lat_F >= lat_plot_south);
lat_u_ind = find(lat_u <= lat_plot_north & lat_u >= lat_plot_south);
lat_v_ind = find(lat_v <= lat_plot_north & lat_v >= lat_plot_south);

%
data.SC_u = SC_u(lat_u_ind,lon_u_ind,p_ind);
data.SC_v = SC_v(lat_v_ind,lon_v_ind,p_ind);
data.div_UV_prime = div_UV_prime;

% plot maps at some depths
close all
fig1 = figure(1);

set(gcf,'units','normalized','outerposition',[0 0 1 1])

rowN = 1; colN = 1;
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (flipud(repmat((1:rowN)',1,colN)))';
gap_w = 0.02; % gap width between subplots
gap_h = 0.06; % gap height between subplots
marg_b = 0.03; % bottom margin
marg_t = 0.03; % top margin
marg_l = 0.02; % left margin
marg_r = 0.005; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length

h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

axes(h_axes_sp(1))

% new colormap routine !
levels = 11;
Reds = othercolor('Reds9', levels);
Blues = othercolor('Blues9', levels);
[Reds(1,:), Blues(1,:)] = deal([1 1 1]);
Blues = flipud(Blues);
BluesReds = [Blues; Reds];
magnif = 1000;

dn = 'div_UV_prime';
dn_units = ' ';

if strcmp(dn_units, 'm^2/s')
    Reds_cont = ...
        [0 0.05 0.1 0.5 1 5 10 20 50 100 200 300]*magnif;
else
    Reds_cont = ...
        [0 0.005 0.01 0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5]*magnif;
end

Blues_cont = -fliplr(Reds_cont);

BluesReds_cont = [Blues_cont, Reds_cont(2:end)];
BluesReds_cont_length = length(BluesReds_cont);

cmap_custom = cmapcust(BluesReds,BluesReds_cont);
colormap(cmap_custom);

h_pcolor = pcolor(...
    lon_F(lon_F_ind)-0.0625, ...
    lat_F(lat_F_ind)+0.0625, ...
    data.(dn)(lat_F_ind,lon_F_ind)*magnif);

set(h_pcolor,'linewidth',0.1)

% shading flat
caxis([Blues_cont(1) Reds_cont(end)]);
freezeColors

hold on
% west
plot([SC_lon_u(1,1) SC_lon_u(1,1)], ...
    [SC_lat_v_north(p_ind,1) SC_lat_v_south(p_ind,1)], ...
    'g','linewidth',0.5,'linestyle','-')
% north
plot(SC_lon_u_repelem(p_ind,:), ...
    SC_lat_v_north_repelem(p_ind,:), ...
    'g','linewidth',0.5,'linestyle','-')
% south
plot(SC_lon_u_repelem(p_ind,:), ...
    SC_lat_v_south_repelem(p_ind,:), ...
    'g','linewidth',0.5,'linestyle','-')
% east
plot([SC_lon_u(1,end) SC_lon_u(1,end)], ...
    [SC_lat_v_north(p_ind,end) SC_lat_v_south(p_ind,end)], ...
    'g','linewidth',0.5,'linestyle','-')

hold on
[lat_u_g_prime_NaN_ind, lon_u_g_prime_NaN_ind] = ...
    find(isnan(u_g_prime(lat_u_ind,lon_u_ind,p_ind)));

lat_u_zoom = lat_u(lat_u_ind);
lon_u_zoom = lon_u(lon_u_ind);

lat_u_g_prime_NaN = lat_u_zoom(lat_u_g_prime_NaN_ind);
lon_u_g_prime_NaN = lon_u_zoom(lon_u_g_prime_NaN_ind);

scatter(lon_u_g_prime_NaN, lat_u_g_prime_NaN, 4, ...
    'o', 'k');

hold on
[lat_v_g_prime_NaN_ind, lon_v_g_prime_NaN_ind] = ...
    find(isnan(v_g_prime(lat_v_ind,lon_v_ind,p_ind)));

lat_v_zoom = lat_v(lat_v_ind);
lon_v_zoom = lon_v(lon_v_ind);

lat_v_g_prime_NaN = lat_v_zoom(lat_v_g_prime_NaN_ind);
lon_v_g_prime_NaN = lon_v_zoom(lon_v_g_prime_NaN_ind);

scatter(lon_v_g_prime_NaN, lat_v_g_prime_NaN, 4, ...
    'o', 'k');

hold on
nn = 1;

% [lon_mg, lat_mg]=meshgrid(lon_u(lon_u_ind),lat_u(lat_u_ind));
% v_g_0 = zeros(size(data.SC_u));
% quivers(...
%     lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
%     data.SC_u(1:nn:end, 1:nn:end), ...
%     v_g_0(1:nn:end, 1:nn:end), ...
%     quiv_S, 1, 'm/s', 'r');

[lon_mg, lat_mg]=meshgrid(lon_v(lon_v_ind),lat_v(lat_v_ind));
u_g_0 = zeros(size(data.SC_v));
quivers(...
    lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
    u_g_0(1:nn:end, 1:nn:end), ...
    data.SC_v(1:nn:end, 1:nn:end), ...
    quiv_S, 1, 'm/s', 'b');

% [SC_u_interp2, SC_v_interp2] = ...
%     deal(NaN(length(lat_u_zoom), length(lon_v_zoom), length(pres_mid)));
% for kk = 1 : length(pres_mid)
%     SC_u_interp2(:,:,kk) = ...
%         interp2(lon_u_zoom, lat_u_zoom, data.SC_u, lon_v_zoom, lat_u_zoom);
%     SC_v_interp2(:,:,kk) = ...
%         interp2(lon_v_zoom, lat_v_zoom, data.SC_v, lon_v_zoom, lat_u_zoom);
% end
% [lon_mg, lat_mg]=meshgrid(lon_v_zoom,lat_u_zoom);
% quivers(...
%     lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
%     SC_u_interp2(1:nn:end, 1:nn:end, p_ind), ...
%     SC_v_interp2(1:nn:end, 1:nn:end, p_ind), ...
%     quiv_S, 3, 'm/s', 'k');

if strcmp(dn_units, 'm^2/s')
    title([dn ' (' dn_units ')'])
else
    title(['Fluxes at ' num2str(pres_mid_now) ' dbar'])
end

cmap = colormap(BluesReds);
BluesReds_linspace = ...
    linspace(Blues_cont(1), Reds_cont(end), BluesReds_cont_length);
% cbar = colorbar;
% cbarrow
% set(cbar, 'YTick',BluesReds_linspace, 'YTickLabel',BluesReds_cont/magnif);

% grid 
font_size = 10;
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')

% export_fig(fig1, ['Dropbox/SACS_work/figures/m17_define_currents_all/v'...
%     num2str(pres_mid_now)], ...
%     '-m3', '-nocrop')

close


%% plot maps of U and V SBC
close all
fig1 = figure(1);

set(gcf,'units','normalized','outerposition',[0 0 1 1])

rowN = 1; colN = 1;
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (flipud(repmat((1:rowN)',1,colN)))';
gap_w = 0.02; % gap width between subplots
gap_h = 0.06; % gap height between subplots
marg_b = 0.03; % bottom margin
marg_t = 0.03; % top margin
marg_l = 0.02; % left margin
marg_r = 0.005; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length

h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

axes(h_axes_sp(1))

colormap(cmap_custom);
h_pcolor = pcolor(...
    lon_F(lon_F_ind)-0.0625, ...
    lat_F(lat_F_ind)+0.0625, ...
    data.(dn)(lat_F_ind,lon_F_ind)*magnif);

set(h_pcolor,'linewidth',0.1)

% shading flat
caxis([Blues_cont(1) Reds_cont(end)]);
freezeColors

hold on
% west
plot([SC_lon_u(1,1) SC_lon_u(1,1)], ...
    [SC_lat_v_north(p_ind,1) SC_lat_v_south(p_ind,1)], ...
    'g','linewidth',0.5,'linestyle','-')
% north
plot(SC_lon_u_repelem(p_ind,:), ...
    SC_lat_v_north_repelem(p_ind,:), ...
    'g','linewidth',0.5,'linestyle','-')
% south
plot(SC_lon_u_repelem(p_ind,:), ...
    SC_lat_v_south_repelem(p_ind,:), ...
    'g','linewidth',0.5,'linestyle','-')
% east
plot([SC_lon_u(1,end) SC_lon_u(1,end)], ...
    [SC_lat_v_north(p_ind,end) SC_lat_v_south(p_ind,end)], ...
    'g','linewidth',0.5,'linestyle','-')

hold on
[lat_u_g_prime_NaN_ind, lon_u_g_prime_NaN_ind] = ...
    find(isnan(u_g_prime(lat_u_ind,lon_u_ind,p_ind)));

lat_u_zoom = lat_u(lat_u_ind);
lon_u_zoom = lon_u(lon_u_ind);

lat_u_g_prime_NaN = lat_u_zoom(lat_u_g_prime_NaN_ind);
lon_u_g_prime_NaN = lon_u_zoom(lon_u_g_prime_NaN_ind);

scatter(lon_u_g_prime_NaN, lat_u_g_prime_NaN, 4, ...
    'o', 'k');

hold on
[lat_v_g_prime_NaN_ind, lon_v_g_prime_NaN_ind] = ...
    find(isnan(v_g_prime(lat_v_ind,lon_v_ind,p_ind)));

lat_v_zoom = lat_v(lat_v_ind);
lon_v_zoom = lon_v(lon_v_ind);

lat_v_g_prime_NaN = lat_v_zoom(lat_v_g_prime_NaN_ind);
lon_v_g_prime_NaN = lon_v_zoom(lon_v_g_prime_NaN_ind);

scatter(lon_v_g_prime_NaN, lat_v_g_prime_NaN, 4, ...
    'o', 'k');

hold on
nn = 1;
quiv_S = 5;
SBC_U_g_times_depth = SC_u.*depth_h_u;
SBC_U_g = nansum(SBC_U_g_times_depth,3);
SBC_U_g_mask = SBC_U_g==0;
SBC_U0 = SBC_U_g + U_ek;
SBC_U0(SBC_U_g_mask) = 0;
SBC_U = SBC_U0(lat_u_ind,lon_u_ind);
[lon_mg, lat_mg]=meshgrid(lon_u(lon_u_ind),lat_u(lat_u_ind));
v_g_0 = zeros(size(SBC_U));
quivers(...
    lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
    SBC_U(1:nn:end, 1:nn:end), ...
    v_g_0(1:nn:end, 1:nn:end), ...
    quiv_S, 1, 'm/s', 'r');

SBC_V_g_times_depth = SC_v.*depth_h_v;
SBC_V_g_raw = nansum(SBC_V_g_times_depth,3);
SBC_V_g_raw_mask = SBC_V_g_raw==0;
SBC_V_raw0 = SBC_V_g_raw + V_ek;
SBC_V_raw0(SBC_V_g_raw_mask) = 0;
SBC_V_raw = SBC_V_raw0(lat_v_ind,lon_v_ind);

SBC_V = NaN(size(SBC_V_raw));
for jj = 1 : length(lon_v_zoom)
    first_now = find(SBC_V_raw(:,jj)~=0,1,'first');
    last_now = find(SBC_V_raw(:,jj)~=0,1,'last');
    SBC_V([first_now,last_now],jj) = SBC_V_raw([first_now,last_now],jj);
end

[lon_mg, lat_mg]=meshgrid(lon_v(lon_v_ind),lat_v(lat_v_ind));
u_g_0 = zeros(size(SBC_V));
quivers(...
    lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
    u_g_0(1:nn:end, 1:nn:end), ...
    SBC_V(1:nn:end, 1:nn:end), ...
    quiv_S, 1, 'm/s', 'b');

% [SC_u_interp2, SC_v_interp2] = ...
%     deal(NaN(length(lat_u_zoom), length(lon_v_zoom), length(pres_mid)));
% for kk = 1 : length(pres_mid)
%     SC_u_interp2(:,:,kk) = ...
%         interp2(lon_u_zoom, lat_u_zoom, data.SC_u, lon_v_zoom, lat_u_zoom);
%     SC_v_interp2(:,:,kk) = ...
%         interp2(lon_v_zoom, lat_v_zoom, data.SC_v, lon_v_zoom, lat_u_zoom);
% end
% [lon_mg, lat_mg]=meshgrid(lon_v_zoom,lat_u_zoom);
% quivers(...
%     lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
%     SC_u_interp2(1:nn:end, 1:nn:end, p_ind), ...
%     SC_v_interp2(1:nn:end, 1:nn:end, p_ind), ...
%     quiv_S, 3, 'm/s', 'k');

if strcmp(dn_units, 'm^2/s')
    title([dn ' (' dn_units ')'])
else
    title(['Fluxes at ' num2str(pres_mid_now) ' dbar'])
end

cmap = colormap(BluesReds);
BluesReds_linspace = ...
    linspace(Blues_cont(1), Reds_cont(end), BluesReds_cont_length);
% cbar = colorbar;
% cbarrow
% set(cbar, 'YTick',BluesReds_linspace, 'YTickLabel',BluesReds_cont/magnif);

% grid 
font_size = 10;
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')

% export_fig(fig1, ['Dropbox/SACS_work/figures/m17_define_currents_all/V_'...
%     num2str(pres_mid_now)], ...
%     '-m3', '-nocrop')
% 
% close


%%
save('cars_out/aus8_currents', 'aus8_currents')
disp('aus8_currents saved in cars_out/aus8_currents !')

