%% fig 1: map of SBC currents
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_currents'])
load([data_path 'SACS_data/KDau_currents'])


%%
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;

U_prime_ztop_to_zmid = aus8_currents.ztop_to_zmid.U_prime.mean;
V_prime_ztop_to_zmid = aus8_currents.ztop_to_zmid.V_prime.mean;
V_prime_ztop_to_zmid_interp2 = interp2(...
    lon_v, lat_v, V_prime_ztop_to_zmid, lon_u, lat_u);

U_g_prime_zmid_to_zbot = aus8_currents.zmid_to_zbot.U_g_prime.mean;
V_g_prime_zmid_to_zbot = aus8_currents.zmid_to_zbot.V_g_prime.mean;
V_g_prime_zmid_to_zbot_interp2 = interp2(...
    lon_v, lat_v, V_g_prime_zmid_to_zbot, lon_u, lat_u);

lon_u_ALLC_repelem = [...
    aus8_currents.lon_u_ALLC(:,1), ...
    repelem(aus8_currents.lon_u_ALLC(:,2:end-1), 1, 2), ...
    aus8_currents.lon_u_ALLC(:,end)];
lat_v_SBC_centr_repelem = ...
    repelem(aus8_currents.lat_v_SBC_centr, 1, 2);
lat_v_SBC_north_repelem = ...
    repelem(aus8_currents.lat_v_SBC_north, 1, 2);
lat_v_SBC_south_repelem = ...
    repelem(aus8_currents.lat_v_SBC_south, 1, 2);
lat_v_DRC_centr_repelem = ...
    repelem(aus8_currents.lat_v_DRC_centr, 1, 2);
lat_v_DRC_north_repelem = ...
    repelem(aus8_currents.lat_v_DRC_north, 1, 2);
lat_v_DRC_south_repelem = ...
    repelem(aus8_currents.lat_v_DRC_south, 1, 2);

z_top = aus8_currents.z_top;
z_mid = aus8_currents.z_mid;
z_bot = aus8_currents.z_bot;


%% 5) plot maps of U and V SBC
screen_ratio = 0.75;
fig_n = 1;
rowcols = [2 1];
rowcols_size = [14 6]/screen_ratio; % cm
margs = [0.6 1.2 0.6 0.6]/screen_ratio; % cm
gaps = [0.4 0.8]/screen_ratio; % cm
plot_cbar_gap = 0.3/screen_ratio;
cbar_x = 0.2/screen_ratio;
cbar_y = rowcols_size(2);

magnif = 10;
cmap1_cont = -[200 100 50 20 10 5 2 0.5 0]*magnif;
cmap2_cont = -fliplr(cmap1_cont);
lvl_cmap1 = length(cmap1_cont)-1;
lvl_cmap2 = length(cmap2_cont)-1;
cmap1 = flipud(othercolor('Blues8', lvl_cmap1));
cmap2 = othercolor('Reds8', lvl_cmap2);
cmap1(end,:) = [1 1 1];
cmap2(1,:) = [1 1 1];
cmaps = [cmap1; cmap2];
cmaps_cont = [cmap1_cont cmap2_cont(2:end)];
cmaps_cont_length = length(cmaps_cont);
cmaps_linspace = linspace(0,1,cmaps_cont_length);
cmaps_y_label = cmaps_cont/magnif;

% lon_min = 115; lon_max = 147; lat_min = -47; lat_max = -32;
lon_min = 108; lon_max = 153; lat_min = -50; lat_max = -30;

x_chc = {aus8_coor.lon_u, aus8_coor.lon_v};
x_ind = [1 1];
y_chc = {aus8_coor.lat_u, aus8_coor.lat_v};

data = {U_prime_ztop_to_zmid*magnif, U_g_prime_zmid_to_zbot*magnif};
v_data = {V_prime_ztop_to_zmid_interp2*magnif, ...
    V_g_prime_zmid_to_zbot_interp2*magnif};

title_chc = {'U''', 'U_{g}'''};
z1_chc = {z_top, z_mid};
z2_chc = {z_mid, z_bot};

for sp = 1 : rowcols(1)*rowcols(2)
    minmax{sp} = [cmaps_cont(1) cmaps_cont(end)];
    cmaps_custom{sp} = cmapcust(cmaps,cmaps_cont);
    
    axis_setup{sp} = [lon_min lon_max lat_min lat_max];
    x{sp} = x_chc{x_ind(sp)};
    y{sp} = y_chc{x_ind(sp)};

end

lett = 'a':'z';
font_size = 8*screen_ratio;
nan_color = [0.7 0.7 0.7];
fig_color = [1 1 1];

close all
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
fig_x = marg_l+colN*x_sp+gap_w*(colN-1)+marg_r;
fig_y = marg_b+rowN*y_sp+gap_h*(rowN-1)+marg_t;

desired_length = 0.05/screen_ratio; %cm
if y_sp > x_sp, long_side = y_sp; else, long_side = x_sp; end
norm_length = desired_length/long_side;
fig_tick_length = [norm_length; 0.01];
if cbar_y > cbar_x, long_side = cbar_y; else, long_side = cbar_x; end
norm_length = desired_length/long_side;
cbar_tick_length = [norm_length; 0.01];

set(fig,'units','centimeters','paperunits','centimeters', ...
    'inverthardcopy','off','color',fig_color,...
    'paperposition',[0 0 fig_x fig_y]*screen_ratio,...
    'position',[0 0 fig_x fig_y]);

for sp = 1 : rowN*colN
    subplot_x = marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1);
    subplot_y = marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1);
    ax = axes('Units','centimeters', ...
        'Position',[subplot_x,subplot_y,x_sp,y_sp]);
    
    % create reference arrow
    lat_ref = -34; lon_ref = 144;
    lat_ref_ind = find(lat_u==lat_ref+1/16); lon_ref_ind = find(lon_u==lon_ref);
    
    data{sp}(lat_ref_ind-15:lat_ref_ind,lon_ref_ind:lon_ref_ind+32) = 0;
    v_data{sp}(lat_ref_ind-15:lat_ref_ind,lon_ref_ind:lon_ref_ind+32) = 0;
    
    if sp == 1
        s = 3;
        ref_magn = 20;
    else
        s = 3;
        ref_magn = 60;
    end
    
    data{sp}(lat_ref_ind-9,lon_ref_ind+6) = ref_magn*magnif;
    
    colormap(ax, cmaps_custom{sp});
    pcolor(x{sp}, y{sp}, data{sp})
    axis(axis_setup{sp})
    shading interp
    caxis([minmax{sp}(1) minmax{sp}(2)]);
    hold on
    
    n = 6; m = 2;
    quiver(lon_u(1:n:end),lat_u(1:m:end),...
        data{sp}(1:m:end,1:n:end),...
        v_data{sp}(1:m:end,1:n:end),...
        s,'k')
    
    text(lon_u(lon_ref_ind+6), lat_u(lat_ref_ind-9+5), ...
        [num2str(ref_magn) ' $m^{2}/s$'], ...
        'fontsize',font_size)
    
    if sp <= 1
        %
        lon_ind = lon_u_ALLC_repelem >= ...
            aus8_currents.nudging.north.west_GAB.lon_u_east & ...
            lon_u_ALLC_repelem <= ...
            aus8_currents.nudging.north.eGs_wBS.lon_u_west;
        h = plot(lon_u_ALLC_repelem(lon_ind), ...
            lat_v_SBC_north_repelem(lon_ind), ...
            'g', 'linewidth', 1);
        %
        lon_ind = lon_u_ALLC_repelem >= ...
            aus8_currents.nudging.north.eGs_wBS.lon_u_east & ...
            lon_u_ALLC_repelem <= ...
            aus8_currents.nudging.north.sBS_wT.lon_u_west;
        h = plot(lon_u_ALLC_repelem(lon_ind), ...
            lat_v_SBC_north_repelem(lon_ind), ...
            'g', 'linewidth', 1);
    end
    %
    h = plot(lon_u_ALLC_repelem, lat_v_SBC_south_repelem, ...
        'r', 'linewidth', 1);
    
    if sp <= 1
        h = plot(lon_u_ALLC_repelem, lat_v_DRC_south_repelem, ...
            '--m', 'linewidth', 1);
    else
        h = plot(lon_u_ALLC_repelem, lat_v_DRC_south_repelem, ...
            ':m', 'linewidth', 1);
    end
    
    
    if sp == 1
        arrow([139 -33], ...
            [lon_u_ALLC_repelem(320) lat_v_SBC_north_repelem(320)], 4)
        text(139, -33, '$y_{CC}$', 'fontsize',font_size)
        arrow([139 -34], ...
            [lon_u_ALLC_repelem(320) lat_v_SBC_south_repelem(320)], 4)
        text(139, -34, '$y_{int}$', 'fontsize',font_size)
        arrow([139 -35], ...
            [lon_u_ALLC_repelem(320) lat_v_DRC_south_repelem(320)], 4)
        text(139, -35, 'upper $y_{OF}$', 'fontsize',font_size)
    end
    
    if sp == 2
        arrow([135 -33.8], ...
            [131.6 -33.8], 4)
        text(135, -33.8, '$z_{int}$ area', 'fontsize',font_size)
        arrow([139 -35], ...
            [lon_u_ALLC_repelem(320) lat_v_DRC_south_repelem(320)], 4)
        text(139, -35, 'lower $y_{OF}$', 'fontsize',font_size)
    end
    
    h_tit = title(['(' lett(sp) ') $' title_chc{sp} ...
        '$ ($m^{2}/s$) integrated from ' ...
    '$z=' num2str(z1_chc{sp}) '$ to $z=' num2str(z2_chc{sp}) '$ $m$'], ...
    'horizontalalignment','left', 'fontsize',font_size);
    h_tit.Position(1) = axis_setup{sp}(1);
    grid
    set(ax,'layer','top','color',nan_color,...
        'fontsize',font_size,'tickdir','out', ...
        'ticklength',fig_tick_length, ...
        'xtick', lon_min:2:lon_max, 'ytick', lat_min:2:lat_max)
    if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
    if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end
    
    if sp == rowN*colN
        ax = axes('visible', 'off');
        colormap(ax, cmaps);
        cbar = colorbar;
        set(cbar,'ytick',cmaps_linspace, ...
            'YAxisLocation','right','YTickLabel',cmaps_y_label,...
            'fontsize',font_size,'ticklength',cbar_tick_length)
        set(cbar,'units','centimeters','position', [...
            (marg_l+x_sp*(cm(sp))+gap_w*(cm(sp)-1)+plot_cbar_gap), ...
            (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
            cbar_x, ...
            cbar_y]);
        pointycbar(cbar)
    end
end

outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
print(fig, ...
    [figures_path mfilename '/' scriptname(1:3) ...
    '_fig' num2str(fig_n) '_'], ...
    '-dpng', '-r300')
close


%%
play(bird_f_path,[1 (get(bird_f_path, 'SampleRate')*3)]);

