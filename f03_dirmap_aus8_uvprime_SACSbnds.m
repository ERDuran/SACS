%% fig 1: map of SBC currents
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_currents'])
load([data_path 'SACS_data/aus8_figures'])
load([data_path 'SACS_data/KDau_fcrt'])


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

U_prime_ztop_to_zmid_K = KDau_fcrt.ztop_to_zmid.fulu.mean;
V_prime_ztop_to_zmid_K = KDau_fcrt.ztop_to_zmid.fulv.mean;
V_prime_ztop_to_zmid_interp2_K = interp2(...
    lon_v, lat_v, V_prime_ztop_to_zmid_K, lon_u, lat_u);

U_g_prime_zmid_to_zbot_K = KDau_fcrt.zmid_to_zbot.fulu.mean;
V_g_prime_zmid_to_zbot_K = KDau_fcrt.zmid_to_zbot.fulv.mean;
V_g_prime_zmid_to_zbot_interp2_K = interp2(...
    lon_v, lat_v, V_g_prime_zmid_to_zbot_K, lon_u, lat_u);

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
rowcols = [4 1];
rowcols_size = [12 4.2]/screen_ratio; % cm
margs = [1.6 0.3 1.9 0.7]/screen_ratio; % cm
gaps = [0.4 0.7]/screen_ratio; % cm

magnif = 10;
cmap1_cont = -[200 100 60 20 5 2 0.5 0]*magnif;
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
lon_min = 110; lon_max = 152; lat_min = -48; lat_max = -31;

x_chc = {aus8_coor.lon_u, aus8_coor.lon_v, aus8_coor.lon_u, aus8_coor.lon_v};
x_ind = [1 1 1 1];
y_chc = {aus8_coor.lat_u, aus8_coor.lat_v, aus8_coor.lat_u, aus8_coor.lat_v};

data = {U_prime_ztop_to_zmid*magnif, U_g_prime_zmid_to_zbot*magnif, ...
    U_prime_ztop_to_zmid_K*magnif, U_g_prime_zmid_to_zbot_K*magnif};
v_data = {V_prime_ztop_to_zmid_interp2*magnif, ...
    V_g_prime_zmid_to_zbot_interp2*magnif, ...
    V_prime_ztop_to_zmid_interp2_K*magnif, ...
    V_g_prime_zmid_to_zbot_interp2_K*magnif};

title_chc = ...
    {'{\boldmath{$V_{up}$}} (arrows) and {$U_{up}$} (shadings)', ...
    '{\boldmath{$V_{low}$}} (arrows) and {$U_{low}$} (shadings)', ...
    '{\boldmath{$V_{up}$}} (arrows) and {$U_{up}$} (shadings)', ...
    '{\boldmath{$V_{low}$}} (arrows) and {$U_{low}$} (shadings)'};
z1_chc = {z_top, z_mid, z_top, z_mid};
z2_chc = {z_mid, z_bot, z_mid, z_bot};

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

plot_cbar_gap = -1.1/screen_ratio;
cbar_x = x_sp;
cbar_y = 0.2/screen_ratio;

desired_length = 0.05/screen_ratio; %cm
if y_sp > x_sp, long_side = y_sp; else, long_side = x_sp; end
norm_length = desired_length/long_side;
fig_tick_length = [norm_length; 0.01];
if cbar_y > cbar_x, long_side = cbar_y; else, long_side = cbar_x; end
norm_length = desired_length/long_side;
cbar_tick_length = norm_length;

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
    elseif sp == 2
        s = 2.75;
        ref_magn = 60;
    elseif sp == 3
        s = 3;
        ref_magn = 20;
    elseif sp == 4
        s = 3.5;
        ref_magn = 60;
    end
    
    data{sp}(lat_ref_ind-10 ,lon_ref_ind+6) = ref_magn*magnif;
    
    colormap(ax, cmaps_custom{sp});
    pcolor(x{sp}, y{sp}, data{sp})
    axis(axis_setup{sp})
    shading interp
    caxis([minmax{sp}(1) minmax{sp}(2)]);
    hold on
    
    for n = 1 : length(aus8_figures.cross.lon(1,:))
        plot(aus8_figures.cross.lon(1:2,n), ...
            aus8_figures.cross.lat(1:2,n), '--k', 'linewidth',1.2)
    end
    
    n = 7; m = 3;
    quiver(lon_u(1:n:end),lat_u(1:m:end),...
        data{sp}(1:m:end,1:n:end),...
        v_data{sp}(1:m:end,1:n:end),...
        s,'k')
    
    text(lon_u(lon_ref_ind+6), lat_u(lat_ref_ind-9+5), ...
        [num2str(ref_magn) ' $m^{2}/s$'], ...
        'fontsize',font_size)
    
    if sp == 1 || sp == 3
        %
        lon_ind = lon_u_ALLC_repelem >= ...
            aus8_currents.nudging.north.west_GAB.lon_u_east & ...
            lon_u_ALLC_repelem <= ...
            aus8_currents.nudging.north.eGs_wBS.lon_u_west;
        h = plot(lon_u_ALLC_repelem(lon_ind), ...
            lat_v_SBC_north_repelem(lon_ind), ...
            'k', 'linewidth', 1.25);
        h = plot(lon_u_ALLC_repelem(lon_ind), ...
            lat_v_SBC_north_repelem(lon_ind), ...
            'g-', 'linewidth', 0.75);
        %
        lon_ind = lon_u_ALLC_repelem >= ...
            aus8_currents.nudging.north.eGs_wBS.lon_u_east & ...
            lon_u_ALLC_repelem <= ...
            aus8_currents.nudging.north.sBS_wT.lon_u_west;
        h = plot(lon_u_ALLC_repelem(lon_ind), ...
            lat_v_SBC_north_repelem(lon_ind), ...
            'k', 'linewidth', 1.25);
        h = plot(lon_u_ALLC_repelem(lon_ind), ...
            lat_v_SBC_north_repelem(lon_ind), ...
            'g-', 'linewidth', 0.75);
    end
    %
    h = plot(lon_u_ALLC_repelem, lat_v_SBC_south_repelem, ...
        'k', 'linewidth', 1.25);
    h = plot(lon_u_ALLC_repelem, lat_v_SBC_south_repelem, ...
        'r-', 'linewidth', 0.75);
    
    if sp == 1
        h = plot(lon_u_ALLC_repelem, lat_v_DRC_south_repelem, ...
            'k', 'linewidth', 1.25);
        h = plot(lon_u_ALLC_repelem, lat_v_DRC_south_repelem, ...
            'm-', 'linewidth', 0.75);
    else
        h = plot(lon_u_ALLC_repelem, lat_v_DRC_south_repelem, ...
            'k', 'linewidth', 1.25);
        h = plot(lon_u_ALLC_repelem, lat_v_DRC_south_repelem, ...
            'm-', 'linewidth', 0.75);
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
        text(139, -35, '$y_{OF}$', 'fontsize',font_size)
    end
    
    if sp == 2
        arrow([135 -33.8], ...
            [131.6 -33.8], 4)
        text(135, -33.8, '$z_{int}$ area', 'fontsize',font_size)
    end
    
    if sp == 3
        small_arr = 3;
        % Cape Leeuwin
        arrow([117,-32], [115,-34.15], small_arr, ...
            'Facecolor', 'k', 'edgecolor', 'k')
        text(117,-32,'Cape Leeuwin', 'fontsize', font_size, ...
            'color', 'k')
        % Albany
        arrow([119,-33.5], [118,-34.9], small_arr, ...
            'Facecolor', 'k', 'edgecolor', 'k')
        text(119,-33.5,'Albany', 'fontsize', font_size, ...
            'color', 'k')
        % Cape Pasley
        arrow([124,-32], [123.4,-33.85], small_arr, ...
            'Facecolor', 'k', 'edgecolor', 'k')
        text(124,-32,'Cape Pasley', 'fontsize', font_size, ...
            'color', 'k')
        % Great Australian Bight
        arrow([132.5,-31.5], [130,-32.3], small_arr, ...
            'Facecolor', 'k', 'edgecolor', 'k')
        text(132.5,-31.5,'Great Australian Bight', 'fontsize', font_size, ...
            'color', 'k')
        % Cape Carnot
        arrow([135.9,-32.9], [135.7,-34.7], small_arr, ...
            'Facecolor', 'k', 'edgecolor', 'k')
        text(135.9,-32.9,sprintf('Cape Carnot'), 'fontsize', font_size, ...
            'color', 'k')
        % Spencer Gulf
%         arrow([142,-31.5], [136.9,-34.2], small_arr, ...
%             'Facecolor', [1 1 1], 'edgecolor', [1 1 1])
%         text(142,-31.5,'Spencer Gulf', 'fontsize', font_size, ...
%             'color', 'k')
        % St. Vincent Gulf
%         arrow([142,-32.5], [138.1,-35], small_arr, ...
%             'Facecolor', [1 1 1], 'edgecolor', [1 1 1])
%         text(142,-32.5,'St. Vincent Gulf', 'fontsize', font_size, ...
%             'color', [1 1 1])
        % Kangaroo Is.
        % arrow([142,-34], [137,-35.85], small_arr, ...
        %     'Facecolor', [1 1 1], 'edgecolor', [1 1 1])
        % text(142,-34,sprintf('Kangaroo\nIs.'), 'fontsize', font_size, ...
        %     'color', [1 1 1])
        % Cape Jaffa
        arrow([141,-35], [139.7,-37], small_arr, ...
            'Facecolor', 'k', 'edgecolor', 'k')
        text(141,-35,sprintf('Cape\nJaffa'), 'fontsize', font_size, ...
            'color', 'k')
        % Portland
        arrow([144.2,-35], [141.3,-38], small_arr, ...
            'Facecolor', 'k', 'edgecolor', 'k')
        text(144.2,-35,'Portland', 'fontsize', font_size, ...
            'color', 'k')
        % Cape Otway
        % arrow([145,-36], [143.5,-38.75], small_arr, ...
        %     'Facecolor', [1 1 1], 'edgecolor', [1 1 1])
        % text(145,-36,'Cape Otway', 'fontsize', font_size, ...
        %     'color', [1 1 1])
        % Bass Strait
        arrow([148.4,-33.5], [145.6,-39.4], small_arr, ...
            'Facecolor', 'k', 'edgecolor', 'k')
        text(148.4,-33.5,sprintf('Bass\nStrait'), 'fontsize', font_size, ...
            'color', 'k')
        % King Is.
        arrow([144.2,-37.5], [144,-39.8], small_arr, ...
            'Facecolor', [0 0 0], 'edgecolor', [0 0 0])
        text(144.2,-37.5,sprintf('King\nIs.'), 'fontsize', font_size, ...
            'color', [0 0 0])
        % Schouten Is.
%         arrow([149.3,-40.8], [148.2 -42.3], small_arr, ...
%             'Facecolor', [0 0 0], 'edgecolor', [0 0 0])
%         text(149.3,-40.8,sprintf('Schouten\nIs.'), 'fontsize', font_size, ...
%             'color', [0 0 0])
        % South East Cape
        arrow([147.9,-36.5], [146.7,-43.62], small_arr, ...
            'Facecolor', [0 0 0], 'edgecolor', [0 0 0])
        text(147.9,-36.5,sprintf('South\nEast\nCape'), 'fontsize', font_size, ...
            'color', [0 0 0])
        % Tasman Rise
        % text(149.5,-47,sprintf('Tasman\nRise'), 'fontsize', font_size, ...
        %     'color', [0 0 0])
        % Tasman Sea
%         arrow([151,-43], [151.9,-43], small_arr, ...
%             'Facecolor', [0 0 0], 'edgecolor', [0 0 0])
%         text(151.5,-43,sprintf('Tasman\nSea'), 'fontsize', font_size, ...
%             'HorizontalAlignment','right','color', [0 0 0])
        % South Australian Basin
        arrow([143.5,-31.5], [129,-38.5], small_arr, ...
            'Facecolor', [0 0 0], 'edgecolor', [0 0 0])
        text(143.5,-31.5,sprintf('South Australian Basin'), ...
            'fontsize', font_size, ...
            'color', [0 0 0])
    end
    
    
    
    h_tit = title(['(' lett(sp) ') ' title_chc{sp} ...
        ' integrated from ' ...
    '$z=' num2str(z1_chc{sp}) '$ to $z=' num2str(z2_chc{sp}) '$ $m$'], ...
    'horizontalalignment','left', 'fontsize',font_size, ...
    'Interpreter','latex');
    h_tit.Position(1) = axis_setup{sp}(1);
    grid
    set(ax,'layer','top','color',nan_color,...
        'fontsize',font_size,'tickdir','out', ...
        'ticklength',fig_tick_length, ...
        'xtick', lon_min:2:lon_max, 'ytick', lat_min:2:lat_max)
    if row_ind(sp) ~= rowN
        set(gca,'xticklabel','')
    else
        xlabel('Longitude ($^{\circ}$E)', 'fontsize', font_size)
    end
    if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end
    
    ylabel('Latitude ($^{\circ}$N)', 'fontsize', font_size)
    
    if sp == rowN*colN
        rectangle('Position',[105 -51 47.75 41.5], ...
            'EdgeColor', [0 0 0], 'linewidth', 1, ...
            'Clipping', 'off')
        text(106.5, -33.5, '\textbf{MOM01-75z}', 'Rotation',90,...
            'fontsize', font_size+2)
        rectangle('Position',[105 -9.0 47.75 39.25], ...
            'EdgeColor', [0 0 0], 'linewidth', 1, ...
            'Clipping', 'off')
        text(106.5, 7, '\textbf{CARS-aus8}', 'Rotation',90,...
            'fontsize', font_size+2)
        
        ax = axes('visible', 'off');
        colormap(ax, cmaps);
        cbar = colorbar('horizontal');
        set(cbar,'xtick',cmaps_linspace, ...
            'XTickLabel',cmaps_y_label,...
            'fontsize',font_size,'ticklength',cbar_tick_length)
        set(cbar,'units','centimeters','position', [...
            marg_l, ...
            marg_b+plot_cbar_gap, ...
            cbar_x, ...
            cbar_y]);
        
        set(get(cbar,'xlabel'), ...
            'String','U ($m^{2}/s$)', ...
            'fontsize',font_size)
        cbar.Label.Interpreter = 'latex';

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
    '_fig' num2str(fig_n) '_2'], ...
    '-dpng', '-r600')
print(fig, ...
    ['~/Duran2017/SACS/10319442jbhpxfsdfvwy/' scriptname(1:3) ...
    '_fig' num2str(fig_n) '_2'], ...
    '-dpng', '-r600')
close


%%
play(bird_f_path,[1 (get(bird_f_path, 'SampleRate')*3)]);

