%% define current boundaries
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/KDau_u_g_prime'])
load([data_path 'SACS_data/KDau_v_g_prime'])
load([data_path 'SACS_data/KDau_U_ek'])
load([data_path 'SACS_data/KDau_V_ek'])
load([data_path 'SACS_data/aus8_currents'])

depth = aus8_coor.depth;
depth_mid = aus8_coor.depth_mid;
depth_thkn = aus8_coor.depth_thkn;
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;
u_bottom_KDau = aus8_coor.u_bottom_KDau;
v_bottom_KDau = aus8_coor.v_bottom_KDau;
u_bottom_KDau(u_bottom_KDau==0) = NaN;
v_bottom_KDau(v_bottom_KDau==0) = NaN;
Months = aus8_coor.Months;
Months{13} = 'mean';


%% 1) Define mid pressure and west, east start/end
%%% BOTTOM PRES
% pressure of the bottom surface currents
z_top = aus8_currents.z_top;
z_mid = aus8_currents.z_mid;
z_bot = aus8_currents.z_bot;
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

lon_u_ALLC = aus8_currents.lon_u_ALLC;
lat_v_SBC_north = aus8_currents.lat_v_SBC_north;
lat_v_SBC_south = aus8_currents.lat_v_SBC_south;
lat_v_DRC_north = aus8_currents.lat_v_DRC_north;
lat_v_DRC_south = aus8_currents.lat_v_DRC_south;


%% 3) Calculate UV_upper and UV_lower
depth_thkn_perm = permute(depth_thkn, [3 2 1]);
depth_thkn_u = repmat(depth_thkn_perm, [length(lat_u), length(lon_u)]);
depth_thkn_v = repmat(depth_thkn_perm, [length(lat_v), length(lon_v)]);

u_g_prime_times_depth_thkn_u = KDau_u_g_prime.mean .* depth_thkn_u;
U_g_prime_ptop_to_pmid = ...
    nansum(u_g_prime_times_depth_thkn_u(:,:,z_mid_above_all_ind), 3);
U_prime_ptop_to_pmid = U_g_prime_ptop_to_pmid + KDau_U_ek.mean;

v_g_prime_times_depth_thkn_v = KDau_v_g_prime.mean .* depth_thkn_v;
V_g_prime_ptop_to_pmid = ...
    nansum(v_g_prime_times_depth_thkn_v(:,:,z_mid_above_all_ind), 3);
V_prime_ptop_to_pmid = V_g_prime_ptop_to_pmid + KDau_V_ek.mean;

U_g_prime_pmid_to_pbot = ...
    nansum(u_g_prime_times_depth_thkn_u(:,:,z_mid_below_all_ind), 3);

V_g_prime_pmid_to_pbot = ...
    nansum(v_g_prime_times_depth_thkn_v(:,:,z_mid_below_all_ind), 3);

U_g_prime_pmid_to_pbot(U_g_prime_pmid_to_pbot==0) = NaN;
V_g_prime_pmid_to_pbot(V_g_prime_pmid_to_pbot==0) = NaN;

KDau_currents.ptop_to_pmid.U_prime.mean = U_prime_ptop_to_pmid;
KDau_currents.ptop_to_pmid.V_prime.mean = V_prime_ptop_to_pmid;
KDau_currents.pmid_to_pbot.U_g_prime.mean = U_g_prime_pmid_to_pbot;
KDau_currents.pmid_to_pbot.V_g_prime.mean = V_g_prime_pmid_to_pbot;


%% 5) make repelem and interp2
% repelem
lon_u_ALLC_repelem = [...
    lon_u_ALLC(:,1), ...
    repelem(lon_u_ALLC(:,2:end-1), 1, 2), ...
    lon_u_ALLC(:,end)];
lat_v_SBC_north_repelem = ...
    repelem(lat_v_SBC_north, 1, 2);
lat_v_SBC_south_repelem = ...
    repelem(lat_v_SBC_south, 1, 2);
lat_v_DRC_north_repelem = ...
    repelem(lat_v_DRC_north, 1, 2);
lat_v_DRC_south_repelem = ...
    repelem(lat_v_DRC_south, 1, 2);


%% 6) plot maps of U and V SBC
close all
fig_n = 1;
rowcols = [2 2];
rowcols_size = [17.5 9]; % cm
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
marg_b = margs(3); % bottom_KDau margin
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
% h = contour(lon_v, lat_u, F_bottom_KDau_depth_now, depth_contours, ...
%     'g', 'linewidth', 0.5);
% % ch = clabel(h, 'manual', 'fontsize', font_size);
% % set(findobj(ch,'String',num2str(SBC_north_p),'String', ...
% %     ['x_{SA}: ' num2str(SBC_north_p)])
% 
% depth_contours = [SBC_south_p SBC_south_p];
% h = contour(lon_v, lat_u, F_bottom_KDau_depth_now, depth_contours, ...
%     'g--', 'linewidth', 0.5);
% % ch = clabel(h, 'manual', 'fontsize', font_size);
% % set(findobj(ch,'String',num2str(SBC_south_p)),'String', ...
% %     ['x_{SBC} : ' num2str(SBC_south_p)])
% 
% depth_contours = [DRC_north_p DRC_north_p];
% h = contour(lon_v, lat_u, F_bottom_KDau_depth_now, depth_contours, ...
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
title(['KDS75 $U''$ ($m^{2/s}$) integrated from ' ...
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
title(['KDS75 $V''$ ($m^{2/s}$) integrated from ' ...
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
title(['KDS75 $U_g''$ ($m^{2/s}$) integrated from ' ...
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
title(['KDS75 $V_g''$ ($m^{2/s}$) integrated from ' ...
    '$z=' num2str(z_mid) '$ to $z=' num2str(z_bot) '$ $m$'])
grid
set(gca,'xtick',lon_min:2:lon_max)
set(gca,'ytick',lat_min:1:lat_max)
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size)
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

%
% outputls = ls(figures_path);
% scriptname = mfilename;
% if ~contains(outputls, scriptname)
%     mkdir(figures_path, scriptname)
% end
% export_fig(fig, ...
%     [figures_path mfilename '/' scriptname(1:3) ...
%     '_fig' num2str(fig_n) '_'], ...
%     '-m3')
% close


%% STUFF
% stuff...
% ) contours
% depth_contours = -700;
% contour(lon_v, lat_v, v_bottom_KDau, [depth_contours depth_contours], ...
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
save([data_path 'SACS_data/KDau_currents'], 'KDau_currents')
disp('KDau_currents DONE')

