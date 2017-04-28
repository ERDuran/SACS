%% define current boundaries
clearvars('-except', '*_path')

load aus8_ZD_method
load topog_mask
pres = (0:10:2000)';
pres_mid = aus8_ZD_method.pres_mid;
depth_h_raw = aus8_ZD_method.depth_thicknesses;
lat_u = aus8_ZD_method.lat_u;
lon_u = aus8_ZD_method.lon_u;
u_g_prime = aus8_ZD_method.u_g_prime;
lat_v = aus8_ZD_method.lat_v;
lon_v = aus8_ZD_method.lon_v;
v_g_prime = aus8_ZD_method.v_g_prime;
lat_F = aus8_ZD_method.lat_F;
lon_F = aus8_ZD_method.lon_F;
div_UV_prime = aus8_ZD_method.div_UV_prime;
U_ek = aus8_ZD_method.U_ek;
V_ek = aus8_ZD_method.V_ek;
F_bottom_depth = round(aus8_ZD_method.F_bottom_depth/10)*10;


%% 1) Define mid pressure and west, east start/end
%%% BOTTOM PRES
% pressure of the bottom surface currents
p_top = 0;
p_mid = 150;
p_bot = 2000;
% pressure above the bottom
p_top_below = p_top+5;
p_mid_above = p_mid-5;
p_mid_below = p_mid+5;
% index vector of the bottom surface currents
p_top_below_ind = find(pres_mid == p_top_below);
p_mid_above_ind = find(pres_mid == p_mid_above);
p_mid_below_ind = find(pres_mid == p_mid_below);
% index vector of currents pressure from surface to interface
p_mid_above_all_ind = pres_mid <= p_mid_above;
p_mid_below_all_ind = pres_mid >= p_mid_below;
%%% BOTTOM PRES


%%% WEST & EAST AND LON UV IND
% longitude of western boundary of the control volume
ALLC_lon_u_west = 115;
% index of longitude of western boundary of the control volume
ALLC_lon_u_west_ind = find(lon_u == ALLC_lon_u_west);
% longitude of eastern boundary of the control volume
ALLC_lon_u_east = 147;
% index of longitude of eastern boundary of the control volume
ALLC_lon_u_east_ind = find(lon_u == ALLC_lon_u_east);
% vector ind
ALLC_lon_u = ALLC_lon_u_west : 1/8 : ALLC_lon_u_east;
ALLC_lon_u_ind = ALLC_lon_u_west_ind : ALLC_lon_u_east_ind;
ALLC_lon_v = ALLC_lon_u_west+1/16 : 1/8 : ALLC_lon_u_east-1/16;
ALLC_lon_v_ind = ALLC_lon_u_west_ind : ALLC_lon_u_east_ind-1;
%%% WEST & EAST LON U


%%% SAVE
aus8_currents.p_top = p_top;
aus8_currents.p_mid = p_mid;
aus8_currents.p_bot = p_bot;
aus8_currents.ALLC_lon_u = ALLC_lon_u;
aus8_currents.ALLC_lon_u_ind = ALLC_lon_u_ind;
aus8_currents.ALLC_lon_v = ALLC_lon_v;
aus8_currents.ALLC_lon_v_ind = ALLC_lon_v_ind;
%%% SAVE


%% 2) Calculate UV_upper and UV_lower
depth_h = permute(depth_h_raw, [3 2 1]);
depth_h_u = repmat(depth_h, [length(lat_u), length(lon_u)]);
depth_h_v = repmat(depth_h, [length(lat_v), length(lon_v)]);

u_g_prime_times_depth_h_u = u_g_prime .* depth_h_u;
U_g_prime_ptop_to_pmid = ...
    nansum(u_g_prime_times_depth_h_u(:,:,p_mid_above_all_ind), 3);
U_prime_ptop_to_pmid = U_g_prime_ptop_to_pmid + U_ek;

v_g_prime_times_depth_h_v = v_g_prime .* depth_h_v;
V_g_prime_ptop_to_pmid = ...
    nansum(v_g_prime_times_depth_h_v(:,:,p_mid_above_all_ind), 3);
V_prime_ptop_to_pmid = V_g_prime_ptop_to_pmid + V_ek;

U_g_prime_pmid_to_pbot = ...
    nansum(u_g_prime_times_depth_h_u(:,:,p_mid_below_all_ind), 3);

V_g_prime_pmid_to_pbot = ...
    nansum(v_g_prime_times_depth_h_v(:,:,p_mid_below_all_ind), 3);

aus8_currents.ptop_to_pmid.U_prime = U_prime_ptop_to_pmid;
aus8_currents.ptop_to_pmid.V_prime = V_prime_ptop_to_pmid;
aus8_currents.pmid_to_pbot.U_g_prime = U_g_prime_pmid_to_pbot;
aus8_currents.pmid_to_pbot.V_g_prime = V_g_prime_pmid_to_pbot;


%% 3) Get all UV within SBC "raw" region
F_bottom_depth(F_bottom_depth == 1990) = 2000;


%%% LAT V BOTTOM EXTENT
% southern extent of the surface current box, taken from the first land
% at the interface between the surface and the subsurface currents
lat_v_south_extent_at_p_mid = 1.0; % degree latitude
% index vector of the southern extent line at the interface between the SC
% SUBSC
lat_v_south_extent_at_p_mid_vec = ...
    1/8 : 1/8 : lat_v_south_extent_at_p_mid;
lat_u_south_extent_at_p_mid_vec = ...
    lat_v_south_extent_at_p_mid_vec(1:end-1);
% length of the southern extent line at the interface between the SC SUBSC
lat_v_south_extent_at_p_mid_length = ...
    length(lat_v_south_extent_at_p_mid_vec);
lat_u_south_extent_at_p_mid_length = ...
    length(lat_u_south_extent_at_p_mid_vec);
% vector ind
lat_v_south_extent_at_p_mid_vec_ind = ...
    1 : 1 : lat_v_south_extent_at_p_mid_length;
lat_u_south_extent_at_p_mid_vec_ind = ...
    1 : 1 : lat_u_south_extent_at_p_mid_length;
%%% LAT V BOTTOM EXTENT


%%% PREP FOR LAT V AT THE BOTTOM
% templates for:
% latitudes of the interface lines
% indices vectors of the latitudes of the interface lines
SBC_raw_U_prime = NaN(size(U_ek));
SBC_raw_V_prime = NaN(size(V_ek));
%%% PREP FOR LAT V AT THE BOTTOM


%%% GET LAT U AND V AT THE BOTTOM
% for all longitudes within the west and east boundaries
for jj = ALLC_lon_u_ind
    % u g prime north
    u_g_prime_at_p_top_below_now = u_g_prime(:,jj,p_top_below_ind);
    u_g_prime_at_p_top_below_now_first_ocean_ind = ...
        find(isnan(u_g_prime_at_p_top_below_now), 1, 'last') + 1;
    
    % v g prime north
    v_g_prime_at_p_top_below_now = v_g_prime(:,jj,p_top_below_ind);
    v_g_prime_at_p_top_below_now_first_ocean_ind = ...
        find(isnan(v_g_prime_at_p_top_below_now), 1, 'last') + 1;
    
    % u g prime south
    u_g_prime_at_p_mid_above_now = u_g_prime(:,jj,p_mid_above_ind);
    u_g_prime_at_p_mid_above_now_first_ocean_ind = ...
        find(isnan(u_g_prime_at_p_mid_above_now), 1, 'last') + 1;
    
    % v g prime south
    v_g_prime_at_p_mid_above_now = v_g_prime(:,jj,p_mid_above_ind);
    v_g_prime_at_p_mid_above_now_first_ocean_ind = ...
        find(isnan(v_g_prime_at_p_mid_above_now), 1, 'last') + 1;
    
    % north
    if v_g_prime_at_p_top_below_now_first_ocean_ind < ...
            u_g_prime_at_p_top_below_now_first_ocean_ind
        
        p_top_below_now_first_ocean_ind = ...
            v_g_prime_at_p_top_below_now_first_ocean_ind;
    else
        p_top_below_now_first_ocean_ind = ...
            u_g_prime_at_p_top_below_now_first_ocean_ind;
    end
    
    % south
    if v_g_prime_at_p_mid_above_now_first_ocean_ind < ...
            u_g_prime_at_p_mid_above_now_first_ocean_ind
        
        p_mid_above_now_first_ocean_ind = ...
            v_g_prime_at_p_mid_above_now_first_ocean_ind;
    else
        p_mid_above_now_first_ocean_ind = ...
            u_g_prime_at_p_mid_above_now_first_ocean_ind;
    end
    lat_u_p_mid_above_now_last_SBC_ind = ...
        p_mid_above_now_first_ocean_ind + ...
        lat_u_south_extent_at_p_mid_vec_ind(end) - 1;
    lat_v_p_mid_above_now_last_SBC_ind = ...
        p_mid_above_now_first_ocean_ind + ...
        lat_v_south_extent_at_p_mid_vec_ind(end) - 1;
        
    
    % 
    SBC_raw_lat_u_now_ind = ...
        (p_top_below_now_first_ocean_ind : ...
        lat_u_p_mid_above_now_last_SBC_ind)';
    SBC_raw_U_prime...
        (p_top_below_now_first_ocean_ind : ...
        lat_u_p_mid_above_now_last_SBC_ind, jj) = ...
        U_prime_ptop_to_pmid(SBC_raw_lat_u_now_ind,jj);
    
    % 
    if jj ~= ALLC_lon_u_ind(end)
        SBC_raw_lat_v_now_ind = ...
            (p_top_below_now_first_ocean_ind : ...
            lat_v_p_mid_above_now_last_SBC_ind)';
        SBC_raw_V_prime...
            (p_top_below_now_first_ocean_ind : ...
            lat_v_p_mid_above_now_last_SBC_ind, jj) = ...
            V_prime_ptop_to_pmid(SBC_raw_lat_v_now_ind,jj);
        
    end
end
%%% GET LAT U AND V AT THE BOTTOM


%%% SAVE
aus8_currents.SBC_raw.U_prime = SBC_raw_U_prime;
aus8_currents.SBC_raw.V_prime = SBC_raw_V_prime;
%%% SAVE


%% 4) make repelem and interp2
% repelem
ALLC_lon_u_repelem = [...
    ALLC_lon_u(:,1), ...
    repelem(ALLC_lon_u(:,2:end-1), 1, 2), ...
    ALLC_lon_u(:,end)];

[ALLC_lat_v_north, ALLC_lat_v_south] = deal(NaN(1, length(ALLC_lon_u)-1));
jj_count = 0;
for jj = 1 : length(lon_v)
    ALLC_lat_v_north_now = ...
        find(isfinite(SBC_raw_V_prime(:,jj)), 1, 'first');
    ALLC_lat_v_south_now = ...
        find(isfinite(SBC_raw_V_prime(:,jj)), 1, 'last');
    
    if ~isempty(ALLC_lat_v_north_now)
        jj_count = jj_count + 1;
        ALLC_lat_v_north(jj_count) = lat_v(ALLC_lat_v_north_now);
        ALLC_lat_v_south(jj_count) = lat_v(ALLC_lat_v_south_now);
    end
end
ALLC_lat_v_north_repelem = ...
    repelem(ALLC_lat_v_north, 1, 2);
ALLC_lat_v_south_repelem = ...
    repelem(ALLC_lat_v_south, 1, 2);

% interp2
U_prime_interp2 = interp2(...
    lon_u, lat_u, U_prime_ptop_to_pmid, lon_v, lat_u);
V_prime_interp2 = interp2(...
    lon_v, lat_v, V_prime_ptop_to_pmid, lon_v, lat_u);
speed_interp2 = sqrt(U_prime_interp2.^2 + V_prime_interp2.^2);

aus8_currents.SBC_raw.ALLC_lon_u_repelem = ALLC_lon_u_repelem;
aus8_currents.SBC_raw.ALLC_lat_v_north_repelem = ALLC_lat_v_north_repelem;
aus8_currents.SBC_raw.ALLC_lat_v_south_repelem = ALLC_lat_v_south_repelem;
% aus8_currents.SBC_raw.U_prime_interp2 = U_prime_interp2;
% aus8_currents.SBC_raw.V_prime_interp2 = V_prime_interp2;
% aus8_currents.SBC_raw.speed_interp2 = speed_interp2;


%% plot maps of U and V SBC
% 1) figure set-up
close all
font_size = 8;
fig1 = figure(1);
set(gcf, 'units', 'normalized', 'outerposition', [0 0 0.8 1]);
rowN = 1; colN = 1;
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = - 0.04; % gap width between subplots
gap_h = 0.04; % gap height between subplots
marg_b = 0.04; % bottom margin
marg_t = 0.04; % top margin
marg_l = 0.02; % left margin
marg_r = 0.02; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length
h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);
colormap(flipud(othercolor('Reds9')))

% 2) all data set-up
sp = 1;
axes(h_axes_sp(sp))
U_prime_ptop_to_pmid_now = U_prime_ptop_to_pmid;
U_prime_ptop_to_pmid_now(U_prime_ptop_to_pmid_now<0) = 0;
V_prime_ptop_to_pmid_now = V_prime_ptop_to_pmid;
V_prime_ptop_to_pmid_now(V_prime_ptop_to_pmid_now>0) = 0;
U_prime_interp2_now = U_prime_interp2;
U_prime_interp2_now(U_prime_interp2_now<0) = 0;
V_prime_interp2_now = V_prime_interp2;
V_prime_interp2_now(V_prime_interp2_now>0) = 0;
speed_interp2_now = sqrt(U_prime_interp2_now.^2 + V_prime_interp2_now.^2);

% 3) asal pcolor set-up
% new colormap routine !
levels = 11;
Reds = othercolor('Reds9', levels);
Reds(1,:) = [1 1 1];
magnif = 100;
Reds_cont = ...
    [0 0.05 0.1 0.5 1 5 10 20 50 75 100 200]*magnif;
Reds_cont_length = length(Reds_cont);
cmap_custom = cmapcust(Reds,Reds_cont);
colormap(cmap_custom);

% 4) plot asal pcolor
% pcolor(...
%     lon_u, ...
%     lat_u, ...
%     U_prime_ptop_to_pmid_now*magnif)
% pcolor(...
%     lon_v, ...
%     lat_v, ...
%     -V_prime_ptop_to_pmid_now*magnif)
pcolor(...
    lon_v, ...
    lat_u, ...
    speed_interp2_now*magnif)
shading interp
hold on
caxis([Reds_cont(1) Reds_cont(end)]);
axis([114 148 -45 -32])
freezeColors

cmap = colormap(Reds);
Reds_linspace = ...
    linspace(Reds_cont(1), Reds_cont(end), Reds_cont_length);
cbar = colorbar;
set(cbar, 'YTick',Reds_linspace, 'YTickLabel',Reds_cont/magnif);

% ) contours
depth_contours = [100 2000];
h = contour(lon_v, lat_u, F_bottom_depth, depth_contours, 'k');

% 6) 
hold on
[lat_u_NaN_ind, lon_u_NaN_ind] = ...
    find(isnan(U_prime_ptop_to_pmid));
lat_u_NaN = lat_u(lat_u_NaN_ind);
lon_u_NaN = lon_u(lon_u_NaN_ind);
scatter(lon_u_NaN, lat_u_NaN, 4, ...
    'o', 'k');

hold on
[lat_v_NaN_ind, lon_v_NaN_ind] = ...
    find(isnan(V_prime_ptop_to_pmid));
lat_v_NaN = lat_v(lat_v_NaN_ind);
lon_v_NaN = lon_v(lon_v_NaN_ind);
scatter(lon_v_NaN, lat_v_NaN, 4, ...
    'o', 'k');

% 7) UV quiver set-up
[lon_mg, lat_mg]=meshgrid(lon_v,lat_u);
quiv_S = 5;
nn = 1;

% 8) plot directional U and V
quiver(...
    lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
    U_prime_interp2(1:nn:end, 1:nn:end), ...
    V_prime_interp2(1:nn:end, 1:nn:end), ...
    quiv_S, 'k');
% reference arrow goes here

% 9) title, grid, background and fonts
title('SBC current')
grid
set(gca,'xtick',115:2:147)
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end


% outputls = ls(figures_path);
% scriptname = mfilename;
% if ~contains(outputls, scriptname)
%     mkdir(figures_path, scriptname)
% end
% export_fig(fig1, [figures_path mfilename '/' scriptname(1:3) '_'], ...
%     '-m4')
% close
%
%
% %% save
% save([cars_out_path 'aus8_ZD_method'], 'aus8_ZD_method')
% disp(['aus8_ZD_method saved in ' ...
%     cars_out_path 'aus8_ZD_method'])

