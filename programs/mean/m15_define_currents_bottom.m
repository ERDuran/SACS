%% define current boundaries
clearvars('-except', '*_path')

load aus8_ZD_method
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


%% define surface currents
%%% BOTTOM PRES
% pressure of the bottom surface currents
SC_bottom_pres = 200;
% pressure above the bottom
SC_bottom_pres_mid = 205;
% index vector of the bottom surface currents
SC_bottom_pres_mid_ind = pres_mid == SC_bottom_pres_mid;
% index vector of currents pressure from surface to interface
SC_bottom_all_pres_mid_ind = pres_mid <= SC_bottom_pres_mid;
%%% BOTTOM PRES


%%% WEST & EAST AND LON UV IND
% longitude of western boundary of the control volume
ALLC_west_lon_u = 115;
% index of longitude of western boundary of the control volume
ALLC_west_lon_u_ind = find(lon_u == ALLC_west_lon_u);
% longitude of eastern boundary of the control volume
ALLC_east_lon_u = 147;
% index of longitude of eastern boundary of the control volume
ALLC_east_lon_u_ind = find(lon_u == ALLC_east_lon_u);
% vector ind
ALLC_west_to_east_lon_u_ind = ALLC_west_lon_u_ind:ALLC_east_lon_u_ind;
aus8_currents.ALLC_west_lon_u = ALLC_west_lon_u;
aus8_currents.ALLC_east_lon_u = ALLC_east_lon_u;
aus8_currents.ALLC_west_to_east_lon_u_ind = ALLC_west_to_east_lon_u_ind;
%%% WEST & EAST LON U


%%% LAT V BOTTOM EXTENT
% southern extent of the surface current box, taken from the first land
% at the interface between the surface and the subsurface currents
SC_v_bottom_south_1 = 0.625; % degree latitude
SC_v_bottom_south_2 = 0.875;
SC_v_bottom_south_3 = 1.125;
SC_v_bottom_south_4 = 1.375;

% length of the southern extent line at the interface between the SC SUBSC
SC_v_bottom_south_length_1 = SC_v_bottom_south_1 / 0.125;
SC_v_bottom_south_length_2 = SC_v_bottom_south_2 / 0.125;
SC_v_bottom_south_length_3 = SC_v_bottom_south_3 / 0.125;
SC_v_bottom_south_length_4 = SC_v_bottom_south_4 / 0.125;
% index vector of the southern extent line at the interface between the SC
% SUBSC
SC_v_bottom_south_vec_1 = 1:SC_v_bottom_south_length_1;
SC_v_bottom_south_vec_2 = 1:SC_v_bottom_south_length_2;
SC_v_bottom_south_vec_3 = 1:SC_v_bottom_south_length_3;
SC_v_bottom_south_vec_4 = 1:SC_v_bottom_south_length_4;
SC_v_bottom_south_length_1 = length(SC_v_bottom_south_vec_1);
SC_v_bottom_south_length_2 = length(SC_v_bottom_south_vec_2);
SC_v_bottom_south_length_3 = length(SC_v_bottom_south_vec_3);
SC_v_bottom_south_length_4 = length(SC_v_bottom_south_vec_4);
%
SC_u_bottom_south_vec_1 = 1:SC_v_bottom_south_length_1-1;
SC_u_bottom_south_vec_2 = 1:SC_v_bottom_south_length_2-1;
SC_u_bottom_south_vec_3 = 1:SC_v_bottom_south_length_3-1;
SC_u_bottom_south_vec_4 = 1:SC_v_bottom_south_length_4-1;
SC_u_bottom_south_length_1 = length(SC_u_bottom_south_vec_1);
SC_u_bottom_south_length_2 = length(SC_u_bottom_south_vec_2);
SC_u_bottom_south_length_3 = length(SC_u_bottom_south_vec_3);
SC_u_bottom_south_length_4 = length(SC_u_bottom_south_vec_4);
%%% LAT V BOTTOM EXTENT


%%% WEST TO EAST LON U AND LON V
% place the u longitudes within the west and east boundaries into the
% template
lon_u_SC_bottom_line = lon_u(ALLC_west_to_east_lon_u_ind);
% lon u grid
SC_lon_u_bottom = ...
    repmat(lon_u_SC_bottom_line, [SC_u_bottom_south_length_4, 1]);
% lon v
lon_v_SC_bottom_line = lon_v(ALLC_west_to_east_lon_u_ind(1:end-1));
SC_lon_v_bottom = ...
    repmat(lon_v_SC_bottom_line,[SC_v_bottom_south_length_4, 1]);
%%% WEST TO EAST LON U AND LON V


%%% PREP FOR LAT V AT THE BOTTOM
% templates for:
% latitudes of the interface lines
% indices vectors of the latitudes of the interface lines
SC_lat_u_bottom = NaN(...
    SC_u_bottom_south_length_4, ...
    length(lon_u_SC_bottom_line));

SC_lat_v_bottom = NaN(...
    SC_v_bottom_south_length_4, ...
    length(lon_v_SC_bottom_line));

SC_lat_v_bottom_south_line = ...
    NaN(1, length(lon_v_SC_bottom_line));

% longitude index counter
jj_count = 0;
%%% PREP FOR LAT V AT THE BOTTOM


%%% GET LAT U AND V AT THE BOTTOM
% for all longitudes within the west and east boundaries
for jj = ALLC_west_to_east_lon_u_ind
    % increase lon counter
    jj_count = jj_count + 1;
    
    % lat v
    SC_bottom_all_lat_u = u_g_prime(:,jj,SC_bottom_pres_mid_ind);
    SC_bottom_lat_u_first_ocean_ind = ...
        find(isnan(SC_bottom_all_lat_u), 1, 'last')+1;
    
    if lon_u(jj) <= 134
        SC_u_bottom_south_vec = SC_u_bottom_south_vec_1;
        SC_v_bottom_south_vec = SC_v_bottom_south_vec_1;
    elseif lon_u(jj) > 134 && lon_u(jj) <= 142
        SC_u_bottom_south_vec = SC_u_bottom_south_vec_2;
        SC_v_bottom_south_vec = SC_v_bottom_south_vec_2;
    elseif lon_u(jj) > 142 && lon_u(jj) <= 142.25
        SC_u_bottom_south_vec = SC_u_bottom_south_vec_3;
        SC_v_bottom_south_vec = SC_v_bottom_south_vec_3;
    elseif lon_u(jj) > 142.25 && lon_u(jj) <= 143.25
        SC_u_bottom_south_vec = SC_u_bottom_south_vec_4;
        SC_v_bottom_south_vec = SC_v_bottom_south_vec_4;
    elseif lon_u(jj) > 143.25 && lon_u(jj) <= 144.75
        SC_u_bottom_south_vec = SC_u_bottom_south_vec_3;
        SC_v_bottom_south_vec = SC_v_bottom_south_vec_3;
    elseif lon_u(jj) > 144.75 && lon_u(jj) <= 145.875
        SC_u_bottom_south_vec = SC_u_bottom_south_vec_2;
        SC_v_bottom_south_vec = SC_v_bottom_south_vec_2;
    elseif lon_u(jj) > 145.875
        SC_u_bottom_south_vec = SC_u_bottom_south_vec_1;
        SC_v_bottom_south_vec = SC_v_bottom_south_vec_1;
    end
    
    SC_bottom_lat_u_ind_now = ...
        (SC_bottom_lat_u_first_ocean_ind+SC_u_bottom_south_vec-1)';
    SC_lat_u_bottom(SC_u_bottom_south_vec,...
        jj_count) = lat_u(SC_bottom_lat_u_ind_now);
    
    if jj ~= ALLC_west_to_east_lon_u_ind(end)
        SC_bottom_lat_v_ind_now = ...
            (SC_bottom_lat_u_first_ocean_ind+SC_v_bottom_south_vec-1)';
        %
        SC_lat_v_bottom(SC_v_bottom_south_vec,...
            jj_count) = lat_v(SC_bottom_lat_v_ind_now);
        
        SC_lat_v_bottom_south_line(jj_count) = ...
            lat_v(SC_bottom_lat_v_ind_now(end));
    end
    
end
%%% GET LAT U AND V AT THE BOTTOM


% SAVE
aus8_currents.SC.bottom.pres = SC_bottom_pres;
aus8_currents.SC.bottom.lat_u = SC_lat_u_bottom;
aus8_currents.SC.bottom.lon_u = SC_lon_u_bottom;
aus8_currents.SC.bottom.lat_v = SC_lat_v_bottom;
aus8_currents.SC.bottom.lon_v = SC_lon_v_bottom;
aus8_currents.SC.bottom.lat_v_south_line = SC_lat_v_bottom_south_line;

% show currents
%
pres_mid_now = SC_bottom_pres_mid;
p_ind = find(pres_mid == pres_mid_now);

%
data.u_g_prime = u_g_prime(:,:,p_ind);
data.v_g_prime = v_g_prime(:,:,p_ind);
data.div_UV_prime = div_UV_prime;

%
lon_plot_west = 114;
lon_plot_east = 149;
lon_F_ind = find(lon_F >= lon_plot_west & lon_F <= lon_plot_east);
lat_plot_north = -31.5;
lat_plot_south = -45.5;
lat_F_ind = find(lat_F <= lat_plot_north & lat_F >= lat_plot_south);

% plot maps at some depths
close all
fig1 = figure(1);

set(gcf,'units','normalized','outerposition',[0 0 2 2])

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

% shading flat
caxis([Blues_cont(1) Reds_cont(end)]);
freezeColors

hold on
lon_u_SC_bottom_grid_plot = SC_lon_u_bottom;
lon_u_SC_bottom_grid_plot(end+1,:) = SC_lon_u_bottom(end,:);
% west
plot(lon_u_SC_bottom_grid_plot(:,1), ...
    SC_lat_v_bottom(:,1), ...
    'g','linewidth',0.5,'linestyle','-')
% north
lon_u_SC_bottom_grid_plot_repelem = [...
    lon_u_SC_bottom_grid_plot(:,1), ...
    repelem(lon_u_SC_bottom_grid_plot(:,2:end-1), 1, 2), ...
    lon_u_SC_bottom_grid_plot(:,end)];
lat_v_SC_bottom_grid_repelem = ...
    repelem(SC_lat_v_bottom, 1, 2);
plot(lon_u_SC_bottom_grid_plot_repelem(1,:), ...
    lat_v_SC_bottom_grid_repelem(1,:), ...
    'g','linewidth',0.5,'linestyle','-')
% south
lat_v_SC_bottom_grid_last_finite_repelem = ...
    repelem(SC_lat_v_bottom_south_line, 1, 2);
plot(lon_u_SC_bottom_grid_plot_repelem(end,:), ...
    lat_v_SC_bottom_grid_last_finite_repelem, ...
    'g','linewidth',0.5,'linestyle','-')
% east
plot(lon_u_SC_bottom_grid_plot(:,end), ...
    SC_lat_v_bottom(:,end-1), ...
    'g','linewidth',0.5,'linestyle','-')

hold on
nn = 1;
quiv_S = 3;
[lon_mg, lat_mg]=meshgrid(lon_u,lat_u);
v_g_0 = zeros(size(u_g_prime));
quiver(...
    lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
    u_g_prime(1:nn:end, 1:nn:end, p_ind), ...
    v_g_0(1:nn:end, 1:nn:end, p_ind), ...
    quiv_S, 'r');

hold on
[lon_mg, lat_mg]=meshgrid(lon_v,lat_v);
u_g_0 = zeros(size(v_g_prime));
quiver(...
    lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
    u_g_0(1:nn:end, 1:nn:end, p_ind), ...
    v_g_prime(1:nn:end, 1:nn:end, p_ind), ...
    quiv_S, 'b');

hold on
[lat_u_g_prime_NaN_ind, lon_u_g_prime_NaN_ind] = ...
    find(isnan(u_g_prime(:,:,p_ind)));

lat_u_g_prime_NaN = lat_u(lat_u_g_prime_NaN_ind);
lon_u_g_prime_NaN = lon_u(lon_u_g_prime_NaN_ind);

scatter(lon_u_g_prime_NaN, lat_u_g_prime_NaN, 4, ...
    'o', 'k');

hold on
[lat_v_g_prime_NaN_ind, lon_v_g_prime_NaN_ind] = ...
    find(isnan(v_g_prime(:,:,p_ind)));

lat_v_g_prime_NaN = lat_v(lat_v_g_prime_NaN_ind);
lon_v_g_prime_NaN = lon_v(lon_v_g_prime_NaN_ind);

scatter(lon_v_g_prime_NaN, lat_v_g_prime_NaN, 4, ...
    'o', 'k');


% plot(SUBC_lon_u_west_south_bound_to_interf, ...
%     SUBC_lat_u_west_south_bound_to_interf, ...
%     'm','linewidth',1.5,'linestyle',':')
% plot(SUBC_lon_u_east_south_bound_to_interf, ...
%     SUBC_lat_u_east_south_bound_to_interf, ...
%     'm','linewidth',1.5,'linestyle',':')


if strcmp(dn_units, 'm^2/s')
    title([dn ' (' dn_units ')'])
else
    title([dn ' (' dn_units ')' ' at ' num2str(pres_mid_now) ' dbar'])
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

% export_fig(fig1, ['figures/m16_define_currents_boundaries/' dn '_' ...
%     num2str(pres_mid_now)], ...
%     '-m3', '-nocrop')

% close


%% save
save([cars_out_path 'aus8_ZD_method'], 'aus8_ZD_method')
disp(['aus8_ZD_method saved in ' ...
    cars_out_path 'aus8_ZD_method'])


