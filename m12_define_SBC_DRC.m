%% define current boundaries
clearvars('-except', '*_path', 'a')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_u_g_prime'])
load([data_path 'SACS_data/aus8_v_g_prime'])
load([data_path 'SACS_data/aus8_U_ek'])
load([data_path 'SACS_data/aus8_V_ek'])

depth = aus8_coor.depth;
depth_mid = aus8_coor.depth_mid;
depth_thkn = aus8_coor.depth_thkn;
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;
u_bottom = aus8_coor.u_bottom;
v_bottom = aus8_coor.v_bottom;
u_bottom(u_bottom==0) = NaN;
v_bottom(v_bottom==0) = NaN;
Months = aus8_coor.Months;
Months{13} = 'mean';


%% 1) Define mid pressure and west, east start/end
%%% BOTTOM PRES
% pressure of the bottom surface currents
z_top = 0;
z_mid = -250;
z_bot = -2000;
% index vector of the bottom surface currents
z_top_below_ind = find(depth==z_top);
z_mid_above_ind = find(depth==z_mid)-1;
z_mid_below_ind = find(depth==z_mid);
% pressure above the bottom
z_top_below = depth_mid(z_top_below_ind);
z_mid_above = depth_mid(z_mid_above_ind);
z_mid_below = depth_mid(z_mid_below_ind);
% index vector of currents pressure from surface to interface
z_mid_above_all_ind = depth_mid >= z_mid_above;
z_mid_below_all_ind = depth_mid <= z_mid_below;
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
lon_u_ALLC = ALLC_lon_u_west : 1/8 : ALLC_lon_u_east;
lon_u_ALLC_ind = ALLC_lon_u_west_ind : ALLC_lon_u_east_ind;
lon_v_ALLC = ALLC_lon_u_west+1/16 : 1/8 : ALLC_lon_u_east-1/16;
lon_v_ALLC_ind = ALLC_lon_u_west_ind : ALLC_lon_u_east_ind-1;
%%% WEST & EAST LON U


%%% SAVE
aus8_currents.z_top = z_top;
aus8_currents.z_mid = z_mid;
aus8_currents.z_bot = z_bot;
aus8_currents.lon_u_ALLC = lon_u_ALLC;
aus8_currents.lon_u_ALLC_ind = lon_u_ALLC_ind;
aus8_currents.lon_v_ALLC = lon_v_ALLC;
aus8_currents.lon_v_ALLC_ind = lon_v_ALLC_ind;
%%% SAVE


%% 2) Define SBC and DRC LATITUDE criteria
% SBC is 
% DRC is 
% IMPORTANT NOTE: if number has four digits, write in the scientific
% form eg. 1000 = 1e+03 or 1700 = 1.7e+03
SBC_centr_zisb = -700; % isobath

DRC_centr_zisb = -2e+03; % isobath

[lat_v_SBC_centr, lat_v_DRC_centr] = deal(NaN(1, length(lon_v_ALLC_ind)));
nc = 0;
for jj = lon_v_ALLC_ind
    nc = nc + 1;
    v_bottom_depth_now = v_bottom(:,jj);
    first_ocean_ind = ...
        find(isnan(v_bottom_depth_now), 1, 'last')+1;
    
    SBC_centr_z_diff = abs(v_bottom_depth_now - SBC_centr_zisb);
%     first_ocean_ind = ...
%         find(isnan(SBC_north_z_diff), 1, 'last')+1;
%     SBC_north_z_diff_first_ocean = ...
%         SBC_north_z_diff(first_ocean_ind:length(lat_v));
%     SBC_north_z_diff_first_pos_ind = ...
%         find(SBC_north_z_diff_first_ocean >= 0, 1, 'first');
%     lat_v_SBC_north(jj_count) = lat_v(...
%         first_ocean_ind-1+SBC_north_z_diff_first_pos_ind-1);
        
    SBC_centr_z_diff_first_ocean = ...
        SBC_centr_z_diff(first_ocean_ind:length(lat_v));
    SBC_centr_z_diff_first_pos_ind = ...
        find(SBC_centr_z_diff_first_ocean == ...
        min(SBC_centr_z_diff_first_ocean), 1, 'first');
    lat_v_SBC_centr(nc) = lat_v(...
        first_ocean_ind-1+SBC_centr_z_diff_first_pos_ind-1);
    
    DRC_centr_z_diff = abs(v_bottom_depth_now - DRC_centr_zisb);
    DRC_centr_z_diff_first_ocean = ...
        DRC_centr_z_diff(first_ocean_ind:length(lat_v));
    DRC_centr_z_diff_first_ocean_min_ind = ...
        find(DRC_centr_z_diff_first_ocean == ...
        min(DRC_centr_z_diff_first_ocean), 1, 'first');
    lat_v_DRC_centr(nc) = lat_v(...
        first_ocean_ind-1+DRC_centr_z_diff_first_ocean_min_ind);
end

% Nudging
% NORTH
SBC_north_dlat1 = +3*1/8; % lat north
SBC_north_dlat2 = +4*1/8; % lat north
lat_v_SBC_north = NaN(1,length(lon_v_ALLC));
WBSEG_lon_u_east = 141.75;
for n = 1 : length(115:1/8:147)-1
    if lon_u_ALLC(n) < WBSEG_lon_u_east
        lat_v_SBC_north(n) = lat_v_SBC_centr(n) + SBC_north_dlat1;
    else
        lat_v_SBC_north(n) = lat_v_SBC_centr(n) + SBC_north_dlat2;
    end
end

% SOUTH
SBC_south_dlat1 = -3*1/8; % lat south
SBC_south_dlat2 = -6*1/8; % lat south
SBC_south_dlat3 = -3*1/8; % lat south
lat_v_SBC_south = NaN(1,length(lon_v_ALLC));
for n = 1 : length(115:1/8:147)-1
    if lon_u_ALLC(n) < WBSEG_lon_u_east
        lat_v_SBC_south(n) = lat_v_SBC_centr(n) + SBC_south_dlat1;
    elseif lon_u_ALLC(n) >= WBSEG_lon_u_east && lon_u_ALLC(n) < 145.75
        lat_v_SBC_south(n) = lat_v_SBC_centr(n) + SBC_south_dlat2;
    else
        lat_v_SBC_south(n) = lat_v_SBC_centr(n) + SBC_south_dlat3;
    end
end

% DRC_north_dlat = SBC_south_dlat
DRC_south_dlat = -1.5; % degree south

lat_v_DRC_north = lat_v_SBC_south;
lat_v_DRC_south = lat_v_DRC_centr + DRC_south_dlat;

%%% WEST GAB AND LON UV IND
% longitude of western boundary of the control volume
lon_u_west = 115;
% index of longitude of western boundary of the control volume
lon_u_west_ind = find(lon_u_ALLC == lon_u_west);
% longitude of eastern boundary of the control volume
lon_u_east = 124;
% index of longitude of eastern boundary of the control volume
lon_u_east_ind = find(lon_u_ALLC == lon_u_east);
% vector ind
lon_u1 = lon_u_west : 1/8 : lon_u_east;
lon_u1_ind = lon_u_west_ind : lon_u_east_ind;
lon_v1 = lon_u_west+1/16 : 1/8 : lon_u_east-1/16;
lon_v1_ind = lon_u_west_ind : lon_u_east_ind-1;
for n = lon_u1_ind(1:end-1)
    v_bottom_depth_now = v_bottom(:,lon_v_ALLC_ind(n));
    u_bottom_depth_now1 = u_bottom(:,lon_u_ALLC_ind(n));
    u_bottom_depth_now2 = u_bottom(:,lon_u_ALLC_ind(n+1));
    if lon_u(lon_u_ALLC_ind(n)) < 116
        v_bottom_depth_now(1:32) = NaN;
        u_bottom_depth_now1(1:32) = NaN;
        u_bottom_depth_now2(1:32) = NaN;
    end
    v_first_ocean_ind = find(isfinite(v_bottom_depth_now), 1, 'first')-1;
    u_first_ocean_ind1 = ...
        find(isfinite(u_bottom_depth_now1), 1, 'first')-1;
    u_first_ocean_ind2 = ...
        find(isfinite(u_bottom_depth_now2), 1, 'first')-1;
    first_ocean = ...
        [v_first_ocean_ind, u_first_ocean_ind1, u_first_ocean_ind2];    
    lat_v_SBC_north(n) = lat_v(min(first_ocean));
end
%%% WEST & EAST LON U

%%% WEST BASS STRAIT, EAST GULFS
% longitude of western boundary of the control volume
lon_u_west = 139.75;
% index of longitude of western boundary of the control volume
lon_u_west_ind = find(lon_u_ALLC == lon_u_west);
% longitude of eastern boundary of the control volume
% WBSEG_lon_u_east = 142;
% index of longitude of eastern boundary of the control volume
lon_u_east_ind = find(lon_u_ALLC == WBSEG_lon_u_east);
% vector ind
lon_u1 = lon_u_west : 1/8 : WBSEG_lon_u_east;
lon_u1_ind = lon_u_west_ind : lon_u_east_ind;
lon_v1 = lon_u_west+1/16 : 1/8 : WBSEG_lon_u_east-1/16;
lon_v1_ind = lon_u_west_ind : lon_u_east_ind-1;
for n = lon_u1_ind(1:end-1)
    v_bottom_depth_now = v_bottom(:,lon_v_ALLC_ind(n));
    u_bottom_depth_now1 = u_bottom(:,lon_u_ALLC_ind(n));
    u_bottom_depth_now2 = u_bottom(:,lon_u_ALLC_ind(n+1));
%     if lon_u(ALLC_lon_u_ind(n)) < 116
%         v_bottom_depth_now(1:32) = NaN;
%         u_bottom_depth_now1(1:32) = NaN;
%         u_bottom_depth_now2(1:32) = NaN;
%     end
    v_first_ocean_ind = find(isnan(v_bottom_depth_now), 1, 'last');
    u_first_ocean_ind1 = ...
        find(isnan(u_bottom_depth_now1), 1, 'last');
    u_first_ocean_ind2 = ...
        find(isnan(u_bottom_depth_now2), 1, 'last');
    first_ocean = ...
        [v_first_ocean_ind, u_first_ocean_ind1, u_first_ocean_ind2];    
    lat_v_SBC_north(n) = lat_v(min(first_ocean));
end
%%% WEST & EAST LON U

%%% WEST TAS AND LON UV IND
% longitude of western boundary of the control volume
lon_u_west = 144.375;
% index of longitude of western boundary of the control volume
lon_u_west_ind = find(lon_u_ALLC == lon_u_west);
% longitude of eastern boundary of the control volume
lon_u_east = 147;
% index of longitude of eastern boundary of the control volume
lon_u_east_ind = find(lon_u_ALLC == lon_u_east);
% vector ind
lon_u1 = lon_u_west : 1/8 : lon_u_east;
lon_u1_ind = lon_u_west_ind : lon_u_east_ind;
lon_v1 = lon_u_west+1/16 : 1/8 : lon_u_east-1/16;
lon_v1_ind = lon_u_west_ind : lon_u_east_ind-1;
for n = lon_u1_ind(1:end-1)
    v_bottom_depth_now = v_bottom(:,lon_v_ALLC_ind(n));
    u_bottom_depth_now1 = u_bottom(:,lon_u_ALLC_ind(n));
    u_bottom_depth_now2 = u_bottom(:,lon_u_ALLC_ind(n+1));
%     if lon_u(ALLC_lon_u_ind(n)) < 116
%         v_bottom_depth_now(1:32) = NaN;
%         u_bottom_depth_now1(1:32) = NaN;
%         u_bottom_depth_now2(1:32) = NaN;
%     end
    v_first_ocean_ind = find(isnan(v_bottom_depth_now), 1, 'last');
    u_first_ocean_ind1 = ...
        find(isnan(u_bottom_depth_now1), 1, 'last');
    u_first_ocean_ind2 = ...
        find(isnan(u_bottom_depth_now2), 1, 'last');
    first_ocean = ...
        [v_first_ocean_ind, u_first_ocean_ind1, u_first_ocean_ind2];    
    lat_v_SBC_north(n) = lat_v(min(first_ocean));
end
%%% WEST & EAST LON U

aus8_currents.lat_v_SBC_north = lat_v_SBC_north;
aus8_currents.lat_v_SBC_south = lat_v_SBC_south;
aus8_currents.lat_v_DRC_north = lat_v_DRC_north;
aus8_currents.lat_v_DRC_south = lat_v_DRC_south;


%% 3) Calculate UV_upper and UV_lower
depth_thkn_perm = permute(depth_thkn, [3 2 1]);
depth_thkn_u = repmat(depth_thkn_perm, [length(lat_u), length(lon_u)]);
depth_thkn_v = repmat(depth_thkn_perm, [length(lat_v), length(lon_v)]);

u_g_prime_times_depth_thkn_u = aus8_u_g_prime.mean .* depth_thkn_u;
U_g_prime_ptop_to_pmid = ...
    nansum(u_g_prime_times_depth_thkn_u(:,:,z_mid_above_all_ind), 3);
U_prime_ptop_to_pmid = U_g_prime_ptop_to_pmid + aus8_U_ek.mean;

v_g_prime_times_depth_thkn_v = aus8_v_g_prime.mean .* depth_thkn_v;
V_g_prime_ptop_to_pmid = ...
    nansum(v_g_prime_times_depth_thkn_v(:,:,z_mid_above_all_ind), 3);
V_prime_ptop_to_pmid = V_g_prime_ptop_to_pmid + aus8_V_ek.mean;

U_g_prime_pmid_to_pbot = ...
    nansum(u_g_prime_times_depth_thkn_u(:,:,z_mid_below_all_ind), 3);

V_g_prime_pmid_to_pbot = ...
    nansum(v_g_prime_times_depth_thkn_v(:,:,z_mid_below_all_ind), 3);

U_g_prime_pmid_to_pbot(U_g_prime_pmid_to_pbot==0) = NaN;
V_g_prime_pmid_to_pbot(V_g_prime_pmid_to_pbot==0) = NaN;

aus8_currents.ptop_to_pmid.U_prime.mean = U_prime_ptop_to_pmid;
aus8_currents.ptop_to_pmid.V_prime.mean = V_prime_ptop_to_pmid;
aus8_currents.pmid_to_pbot.U_g_prime.mean = U_g_prime_pmid_to_pbot;
aus8_currents.pmid_to_pbot.V_g_prime.mean = V_g_prime_pmid_to_pbot;


%% 4) make repelem and interp2
% repelem
lon_u_ALLC_repelem = [...
    lon_u_ALLC(:,1), ...
    repelem(lon_u_ALLC(:,2:end-1), 1, 2), ...
    lon_u_ALLC(:,end)];
lat_v_SBC_centr_repelem = ...
    repelem(lat_v_SBC_centr, 1, 2);
lat_v_SBC_north_repelem = ...
    repelem(lat_v_SBC_north, 1, 2);
lat_v_SBC_south_repelem = ...
    repelem(lat_v_SBC_south, 1, 2);
lat_v_DRC_centr_repelem = ...
    repelem(lat_v_DRC_centr, 1, 2);
lat_v_DRC_north_repelem = ...
    repelem(lat_v_DRC_north, 1, 2);
lat_v_DRC_south_repelem = ...
    repelem(lat_v_DRC_south, 1, 2);


%% 5) plot maps of U and V SBC
close all
fig_n = 1;
rowcols = [2 2];
rowcols_size = [17 9]; % cm
margs = [1 1 1 1]; % cm
gaps = [1 1]; % cm
cmaps_levels = 12;

font_size = 9;
fig_color = [0.7 0.7 0.7];

fig = figure(fig_n);
rowN = rowcols(1); colN = rowcols(2);
[rm, cm] = meshgrid(rowN:-1:1, 1:colN);
x_sp = rowcols_size(1); % x subplot length
y_sp = rowcols_size(2); % y subplot length
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = gaps(1); % gap width between subplots
gap_h = gaps(2); % gap height between subplots
marg_b = margs(3); % bottom margin
marg_t = margs(4); % top margin
marg_l = margs(1); % left margin
marg_r = margs(2); % right margin
set(fig,'units','centimeters',...
    'position',[0 0 ...
    (marg_l+colN*x_sp+gap_w*(colN-1)+marg_r) ...
    (marg_b+rowN*y_sp+gap_h*(rowN-1)+marg_t)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) all data set-up
sp = 1;
axes('Units','centimeters', ...
            'Position',[...
            (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)), ...
            (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
            x_sp, ...
            y_sp])

% 3) asal pcolor set-up
% new colormap routine !
levels = 11;
cmap1 = flipud(othercolor('Blues9', levels));
cmap1(end,:) = [1 1 1];
cmap2 = othercolor('Reds9', levels);
cmap2(1,:) = [1 1 1];
cmaps = [cmap1; cmap2];
magnif = 100;
cmap2_cont = ...
    [0 0.5 1 2 5 7.5 10 20 40 60 80 100]*magnif;
cmap1_cont = -fliplr(cmap2_cont);
cmaps_cont = [cmap1_cont cmap2_cont(2:end)];
cmaps_cont_length = length(cmaps_cont);
cmaps_custom = cmapcust(cmaps,cmaps_cont);
colormap(cmaps_custom);

% 4) plot asal pcolor
pcolor(...
    lon_u, ...
    lat_u, ...
    U_prime_ptop_to_pmid*magnif)
shading interp
hold on
caxis([cmaps_cont(1) cmaps_cont(end)]);
lon_min = 108; lon_max = 152; lat_min = -50; lat_max = -31.5;
axis([lon_min lon_max lat_min lat_max])
freezeColors

cmap = colormap(cmaps);
cmaps_linspace = ...
    linspace(cmaps_cont(1), cmaps_cont(end), cmaps_cont_length);
cbar = colorbar;
set(cbar, 'YTick',cmaps_linspace, 'YTickLabel',cmaps_cont/magnif);

% h = plot(lon_u_ALLC_repelem, lat_v_SBC_centr_repelem, ...
%     'k--', 'linewidth', 0.5);
h = plot(lon_u_ALLC_repelem, lat_v_SBC_north_repelem, ...
    'g', 'linewidth', 2);
h = plot(lon_u_ALLC_repelem, lat_v_SBC_south_repelem, ...
    'r', 'linewidth', 2);
% h = plot(lon_u_ALLC_repelem, lat_v_DRC_centr_repelem, ...
%     'k-', 'linewidth', 0.5);
h = plot(lon_u_ALLC_repelem, lat_v_DRC_north_repelem, ...
    'r', 'linewidth', 2);
h = plot(lon_u_ALLC_repelem, lat_v_DRC_south_repelem, ...
    'k', 'linewidth', 2);


arrow([138 -33], ...
    [lon_u_ALLC_repelem(320) lat_v_SBC_north_repelem(320)])
text(138, -33, 'SBC $V_{t}nc$')
arrow([138 -34], ...
    [lon_u_ALLC_repelem(320) lat_v_SBC_south_repelem(320)])
text(138, -34, 'SBC $V_{t}sc$ and DRC $V_{t}nc$')
arrow([138 -35], ...
    [lon_u_ALLC_repelem(320) lat_v_DRC_south_repelem(320)])
text(138, -35, 'DRC $V_{t}sc$ up')


% ch = clabel(h, 'manual', 'fontsize', font_size);
% set(findobj(ch,'String',num2str(DRC_south_p)),'String', ...
%     ['x_{SBC} : ' num2str(SBC_south_p)])
% depth_contours = [SBC_north_p SBC_north_p];
% h = contour(lon_v, lat_u, F_bottom_depth_now, depth_contours, ...
%     'g', 'linewidth', 0.5);
% % ch = clabel(h, 'manual', 'fontsize', font_size);
% % set(findobj(ch,'String',num2str(SBC_north_p),'String', ...
% %     ['x_{SA}: ' num2str(SBC_north_p)])
% 
% depth_contours = [SBC_south_p SBC_south_p];
% h = contour(lon_v, lat_u, F_bottom_depth_now, depth_contours, ...
%     'g--', 'linewidth', 0.5);
% % ch = clabel(h, 'manual', 'fontsize', font_size);
% % set(findobj(ch,'String',num2str(SBC_south_p)),'String', ...
% %     ['x_{SBC} : ' num2str(SBC_south_p)])
% 
% depth_contours = [DRC_north_p DRC_north_p];
% h = contour(lon_v, lat_u, F_bottom_depth_now, depth_contours, ...
%     'm--', 'linewidth', 0.5);


% 6) 
% hold on
% [lat_u_NaN_ind, lon_u_NaN_ind] = ...
%     find(isnan(U_prime_ptop_to_pmid));
% lat_u_NaN = lat_u(lat_u_NaN_ind);
% lon_u_NaN = lon_u(lon_u_NaN_ind);
% scatter(lon_u_NaN, lat_u_NaN, 4, ...
%     'o', 'w');
% % 
% hold on
% [lat_v_NaN_ind, lon_v_NaN_ind] = ...
%     find(isnan(V_prime_ptop_to_pmid));
% lat_v_NaN = lat_v(lat_v_NaN_ind);
% lon_v_NaN = lon_v(lon_v_NaN_ind);
% scatter(lon_v_NaN, lat_v_NaN, 4, ...
%     'o', 'k');

% 9) title, grid, background and fonts
title(['aus8 $U''$ ($m^{2/s}$) integrated from ' ...
    '$z=0$ to $z=' num2str(z_mid) '$ $m$'])
grid
set(gca,'xtick',lon_min:2:lon_max)
set(gca,'ytick',lat_min:1:lat_max)
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size)
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) all data set-up
sp = 2;
axes('Units','centimeters', ...
            'Position',[...
            (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)), ...
            (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
            x_sp, ...
            y_sp])
% 3) asal pcolor set-up
% new colormap routine !
levels = 11;
cmap1 = flipud(othercolor('Blues9', levels));
cmap1(end,:) = [1 1 1];
cmap2 = othercolor('Reds9', levels);
cmap2(1,:) = [1 1 1];
cmaps = [cmap1; cmap2];
magnif = 100;
cmap2_cont = ...
    [0 0.5 1 2 5 7.5 10 20 40 60 80 100]*magnif;
cmap1_cont = -fliplr(cmap2_cont);
cmaps_cont = [cmap1_cont cmap2_cont(2:end)];
cmaps_cont_length = length(cmaps_cont);
cmaps_custom = cmapcust(cmaps,cmaps_cont);
colormap(cmaps_custom);
% 4) plot asal pcolor
pcolor(...
    lon_v, ...
    lat_v, ...
    V_prime_ptop_to_pmid*magnif)
shading interp
hold on
caxis([cmaps_cont(1) cmaps_cont(end)]);
axis([lon_min lon_max lat_min lat_max])
freezeColors
cmap = colormap(cmaps);
cmaps_linspace = ...
    linspace(cmaps_cont(1), cmaps_cont(end), cmaps_cont_length);
cbar = colorbar;
set(cbar, 'YTick',cmaps_linspace, 'YTickLabel',cmaps_cont/magnif);
% h = plot(lon_u_ALLC_repelem, lat_v_SBC_centr_repelem, ...
%     'k--', 'linewidth', 0.5);
h = plot(lon_u_ALLC_repelem, lat_v_SBC_north_repelem, ...
    'g', 'linewidth', 2);
h = plot(lon_u_ALLC_repelem, lat_v_SBC_south_repelem, ...
    'r', 'linewidth', 2);
% h = plot(lon_u_ALLC_repelem, lat_v_DRC_centr_repelem, ...
%     'k-', 'linewidth', 0.5);
h = plot(lon_u_ALLC_repelem, lat_v_DRC_north_repelem, ...
    'r', 'linewidth', 2);
h = plot(lon_u_ALLC_repelem, lat_v_DRC_south_repelem, ...
    'k', 'linewidth', 2);
% 9) title, grid, background and fonts
title(['aus8 $V''$ ($m^{2/s}$) integrated from ' ...
    '$z=' num2str(z_top) '$ to $z=' num2str(z_mid) '$ $m$'])
grid
set(gca,'xtick',lon_min:2:lon_max)
set(gca,'ytick',lat_min:1:lat_max)
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size)
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) all data set-up
sp = 3;
axes('Units','centimeters', ...
            'Position',[...
            (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)), ...
            (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
            x_sp, ...
            y_sp])
% 3) asal pcolor set-up
% new colormap routine !
levels = 11;
cmap1 = flipud(othercolor('Blues9', levels));
cmap1(end,:) = [1 1 1];
cmap2 = othercolor('Reds9', levels);
cmap2(1,:) = [1 1 1];
cmaps = [cmap1; cmap2];
magnif = 100;
cmap2_cont = ...
    [0 0.5 1 2 5 7.5 10 20 40 60 80 100]*magnif;
cmap1_cont = -fliplr(cmap2_cont);
cmaps_cont = [cmap1_cont cmap2_cont(2:end)];
cmaps_cont_length = length(cmaps_cont);
cmaps_custom = cmapcust(cmaps,cmaps_cont);
colormap(cmaps_custom);
% 4) plot asal pcolor
pcolor(...
    lon_u, ...
    lat_u, ...
    U_g_prime_pmid_to_pbot*magnif)
shading interp
hold on
caxis([cmaps_cont(1) cmaps_cont(end)]);
axis([lon_min lon_max lat_min lat_max])
freezeColors
cmap = colormap(cmaps);
cmaps_linspace = ...
    linspace(cmaps_cont(1), cmaps_cont(end), cmaps_cont_length);
cbar = colorbar;
set(cbar, 'YTick',cmaps_linspace, 'YTickLabel',cmaps_cont/magnif);
% h = plot(lon_u_ALLC_repelem, lat_v_SBC_centr_repelem, ...
%     'k--', 'linewidth', 0.5);
h = plot(lon_u_ALLC_repelem, lat_v_SBC_north_repelem, ...
    'g', 'linewidth', 2);
h = plot(lon_u_ALLC_repelem, lat_v_SBC_south_repelem, ...
    'r', 'linewidth', 2);
% h = plot(lon_u_ALLC_repelem, lat_v_DRC_centr_repelem, ...
%     'k-', 'linewidth', 0.5);
h = plot(lon_u_ALLC_repelem, lat_v_DRC_north_repelem, ...
    'r', 'linewidth', 2);
h = plot(lon_u_ALLC_repelem, lat_v_DRC_south_repelem, ...
    'k', 'linewidth', 2);
arrow([138 -33], ...
    [lon_u_ALLC_repelem(320) lat_v_SBC_north_repelem(320)])
text(138, -33, 'DRC $V_{t}nc$')
arrow([138 -35], ...
    [lon_u_ALLC_repelem(320) lat_v_DRC_south_repelem(320)])
text(138, -35, 'DRC $V_{t}sc$ dw')
% 9) title, grid, background and fonts
title(['aus8 $U_g''$ ($m^{2/s}$) integrated from ' ...
    '$z=' num2str(z_mid) '$ to $z=' num2str(z_bot) '$ $m$'])
grid
set(gca,'xtick',lon_min:2:lon_max)
set(gca,'ytick',lat_min:1:lat_max)
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size)
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) all data set-up
sp = 4;
axes('Units','centimeters', ...
            'Position',[...
            (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)), ...
            (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
            x_sp, ...
            y_sp])
% 3) asal pcolor set-up
% new colormap routine !
levels = 11;
cmap1 = flipud(othercolor('Blues9', levels));
cmap1(end,:) = [1 1 1];
cmap2 = othercolor('Reds9', levels);
cmap2(1,:) = [1 1 1];
cmaps = [cmap1; cmap2];
magnif = 100;
cmap2_cont = ...
    [0 0.5 1 2 5 7.5 10 20 40 60 80 100]*magnif;
cmap1_cont = -fliplr(cmap2_cont);
cmaps_cont = [cmap1_cont cmap2_cont(2:end)];
cmaps_cont_length = length(cmaps_cont);
cmaps_custom = cmapcust(cmaps,cmaps_cont);
colormap(cmaps_custom);
% 4) plot asal pcolor
pcolor(...
    lon_v, ...
    lat_v, ...
    V_g_prime_pmid_to_pbot*magnif)
shading interp
hold on
caxis([cmaps_cont(1) cmaps_cont(end)]);
axis([lon_min lon_max lat_min lat_max])
freezeColors
cmap = colormap(cmaps);
cmaps_linspace = ...
    linspace(cmaps_cont(1), cmaps_cont(end), cmaps_cont_length);
cbar = colorbar;
set(cbar, 'YTick',cmaps_linspace, 'YTickLabel',cmaps_cont/magnif);
% h = plot(lon_u_ALLC_repelem, lat_v_SBC_centr_repelem, ...
%     'k--', 'linewidth', 0.5);
h = plot(lon_u_ALLC_repelem, lat_v_SBC_north_repelem, ...
    'g', 'linewidth', 2);
h = plot(lon_u_ALLC_repelem, lat_v_SBC_south_repelem, ...
    'r', 'linewidth', 2);
% h = plot(lon_u_ALLC_repelem, lat_v_DRC_centr_repelem, ...
%     'k-', 'linewidth', 0.5);
h = plot(lon_u_ALLC_repelem, lat_v_DRC_north_repelem, ...
    'r', 'linewidth', 2);
h = plot(lon_u_ALLC_repelem, lat_v_DRC_south_repelem, ...
    'k', 'linewidth', 2);
% 9) title, grid, background and fonts
title(['aus8 $V_g''$ ($m^{2/s}$) integrated from ' ...
    '$z=' num2str(z_mid) '$ to $z=' num2str(z_bot) '$ $m$'])
grid
set(gca,'xtick',lon_min:2:lon_max)
set(gca,'ytick',lat_min:1:lat_max)
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size)
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

%
outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
export_fig(fig, ...
    [figures_path mfilename '/' scriptname(1:3) ...
    '_fig' num2str(fig_n) '_'], ...
    '-m3')
close


%% STUFF
% stuff...
% ) contours
% depth_contours = -700;
% contour(lon_v, lat_v, v_bottom, [depth_contours depth_contours], ...
%     'k', 'linewidth', 1);
% 
% cn = c;
% cs = c;
% cn(2,2:end) = c(2,2:end)+0.4;
% cs(2,2:end) = c(2,2:end)-0.4;
% % cs(1,2:end) = c(1,2:end)-0.4;
% 
% plot(c(1,2:end), c(2,2:end), 'g', 'linewidth', 1)
% plot(cn(1,2:end), cn(2,2:end), 'g', 'linewidth', 3)
% plot(cs(1,2:end), cs(2,2:end), 'g', 'linewidth', 3)

% 7) UV quiver set-up
% [lon_mg, lat_mg]=meshgrid(lon_v,lat_u);
% qscale = 2/magnif  ; % scaling factor for all vectors
% nn = 1;
% 8) plot directional U and V
% h = quiver(...
%     lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
%     U_prime_interp2(1:nn:end, 1:nn:end), ...
%     V_prime_interp2(1:nn:end, 1:nn:end), ...
%     0, 'k');
% hU = get(h,'UData') ;
% hV = get(h,'VData') ;
% set(h,'UData',qscale*hU,'VData',qscale*hV)
% reference vector at the end

% % reference vector
% ref_vec = 30;
% n_lon = 2;
% lon_length = (lon_max - lon_min);
% lon_x_ratio = x_sp / lon_length;
% lon_x_ratio = lon_x_ratio*n_lon;
% n_lat = 1;
% lat_length = (lat_max - lat_min);
% lat_y_ratio = y_sp / lat_length;
% lat_y_ratio = lat_y_ratio*n_lat;
% axes('Position', ...
%     [marg_l+x_sp-lon_x_ratio-0.06 ...
%     marg_b+y_sp-lat_y_ratio ...
%     lon_x_ratio lat_y_ratio], ...
%     'layer', 'top');
% ref_x = 0:1/8:1*n_lon;
% ref_y = 0:1/8:1*n_lat;
% [ref_x_mg, ref_y_mg] = meshgrid(ref_x, ref_y);
% ref_U = zeros(size(ref_x_mg));
% ref_U(4*n_lat,4*n_lon) = ref_vec;
% ref_V = zeros(size(ref_x_mg));
% 
% h = quiver(...
%     ref_x, ref_y, ...
%     ref_U, ...
%     ref_V, ...
%     0, 'k');
% hU = get(h,'UData') ;
% hV = get(h,'VData') ;
% set(h,'UData',qscale*hU,'VData',qscale*hV)
% set(gca, 'xtick', '')
% set(gca, 'ytick', '')
% text(ref_x(4*n_lon), ref_y(4*n_lat+2), [num2str(ref_vec) ' m^2/s'], ...
%     'fontsize',font_size)


%% 7) save
save([data_path 'SACS_data/aus8_currents'], 'aus8_currents')
disp('aus8_currents DONE')

