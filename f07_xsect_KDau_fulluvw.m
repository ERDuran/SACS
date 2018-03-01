%% fig 2: meridional cross sections of zonal velocity differences
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_currents'])
load([data_path 'SACS_data/aus8_figures'])
load([data_path 'SACS_data/KDau_fulu'])
load([data_path 'SACS_data/KDau_fulv'])
load([data_path 'SACS_data/KDau_fulw'])


%% plot map
depth_mid_for_plot = aus8_coor.depth_mid;
depth_mid_for_plot(end) = -2000;
lon_u_SBC_north_OK = ...
    [aus8_currents.nudging.north.west_GAB.lon_u_east : 1/8 : ...
    aus8_currents.nudging.north.eGs_wBS.lon_u_west, ...
    aus8_currents.nudging.north.eGs_wBS.lon_u_east : 1/8 : ...
    aus8_currents.nudging.north.sBS_wT.lon_u_west];

screen_ratio = 0.75;
fig_n = 1;
rowcols = [2 4];
rowcols_size = [3.3 4]/screen_ratio; % cm
margs = [1.1 1.5 0.8 0.6]/screen_ratio; % cm
gaps = [0.2 1]/screen_ratio; % cm
plot_cbar_gap = 0.2/screen_ratio;
cbar_x = 0.2/screen_ratio;
cbar_y = rowcols_size(2);

magnif = 100;
cmap1_cont = -[14 12 10 8 6 4 2 0.5 0];
cmap2_cont = [0 0.5 2 4 6 8 10 12 14];
lvl_cmap1 = length(cmap1_cont)-1;
lvl_cmap2 = length(cmap2_cont)-1;
cmap1 = flipud(othercolor('Blues6', lvl_cmap1));
cmap2 = othercolor('Reds6', lvl_cmap2);
cmap1(end,:) = [1 1 1];
cmap2(1,:) = [1 1 1];
cmaps = [cmap1; cmap2];
cmaps_cont = [cmap1_cont cmap2_cont(2:end)];
cmaps_cont_length = length(cmaps_cont);
cmaps_linspace = linspace(0,1,cmaps_cont_length);
cmaps_y_label = cmaps_cont;

for sp = 1 : rowcols(1)*rowcols(2)
    minmax{sp} = [cmaps_cont(1) cmaps_cont(end)];
    cmaps_custom{sp} = cmapcust(cmaps,cmaps_cont);
    
    axis_setup{sp} = [...
        aus8_figures.cross.lat(2,sp), ...
        aus8_figures.cross.lat(1,sp), ...
        depth_mid_for_plot(end), ...
        0];
    x{sp} = aus8_coor.lat_u;
    y{sp} = depth_mid_for_plot;
    data{sp} = squeeze(KDau_fulu.mean...
        (:,aus8_coor.lon_u == aus8_figures.cross.lon(1,sp),:))'*magnif;
    data_contour{sp} = squeeze(KDau_fulv.mean...
        (:,aus8_coor.lon_v == aus8_figures.cross.lon(1,sp)+1/16,:))'...
        *magnif;
    x_contour{sp} = aus8_coor.lat_v;
    data_contour_2{sp} = squeeze(KDau_fulw.mean...
        (:,aus8_coor.lon_v == aus8_figures.cross.lon(1,sp)+1/16,:))'...
        *magnif*10000;
    x_contour_2{sp} = aus8_coor.lat_u;
    y_contour_2{sp} = aus8_coor.depth;
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

%
depth_plot = [-20, -100:-100:-2000];
depth_plot_idx = NaN(length(depth_plot),1);
for m = 1 : length(depth_plot)
    depth_plot_idx(m) = ...
        find_nearest(depth_mid_for_plot,depth_plot(m));
end
depth_plot_idx(18) = 61;

for sp = 1 : rowN*colN
    subplot_x = marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1);
    subplot_y = marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1);
    ax = axes('Units','centimeters', ...
        'Position',[subplot_x,subplot_y,x_sp,y_sp]);
    
    colormap(ax, cmaps_custom{sp});
    pcolor(x{sp}, y{sp}, data{sp})
    axis(axis_setup{sp})
    shading interp
    caxis([minmax{sp}(1) minmax{sp}(2)]);
    hold on
    
    [cn, hn] = contour(x_contour{sp},y{sp},data_contour{sp}, ...
        -20 : 2 : -2,'--k');
    clabel(cn, hn,'fontsize',font_size)
    
    [cn, hn] = contour(x_contour{sp},y{sp},data_contour{sp}, ...
        2 : 2 : 20,'-k');
    clabel(cn, hn,'fontsize',font_size)
    
    [cn, hn] = contour(x_contour{sp},y{sp},data_contour{sp}, ...
        [0 0],'-k', 'linewidth', 0.1);
    clabel(cn, hn,'fontsize',font_size)
    
    
    cmap = othercolor('PRGn8',8);
    contour(x_contour_2{sp},y_contour_2{sp}, ...
        data_contour_2{sp}, [-40 -40], 'linecolor', cmap(1,:), ...
        'linewidth', 0.9);
    contour(x_contour_2{sp},y_contour_2{sp}, ...
        data_contour_2{sp}, [-30 -30], 'linecolor', cmap(2,:), ...
        'linewidth', 0.9);
    contour(x_contour_2{sp},y_contour_2{sp}, ...
        data_contour_2{sp}, [-20 -20], 'linecolor', cmap(3,:), ...
        'linewidth', 0.9);
    contour(x_contour_2{sp},y_contour_2{sp}, ...
        data_contour_2{sp}, [-10 -10], 'linecolor', cmap(4,:), ...
        'linewidth', 0.9);
    contour(x_contour_2{sp},y_contour_2{sp}, ...
        data_contour_2{sp}, [10 10], 'linecolor', cmap(5,:), ...
        'linewidth', 0.9);
    contour(x_contour_2{sp},y_contour_2{sp}, ...
        data_contour_2{sp}, [20 20], 'linecolor', cmap(6,:), ...
        'linewidth', 0.9);
    contour(x_contour_2{sp},y_contour_2{sp}, ...
        data_contour_2{sp}, [30 30], 'linecolor', cmap(7,:), ...
        'linewidth', 0.9);
    contour(x_contour_2{sp},y_contour_2{sp}, ...
        data_contour_2{sp}, [40 40], 'linecolor', cmap(8,:), ...
        'linewidth', 0.9);
    
    line_width = 1;
    if find(ismember(lon_u_SBC_north_OK,aus8_figures.cross.lon(1,sp)))
        plot([...
            aus8_currents.lat_v_SBC_north(...
            aus8_currents.lon_u_ALLC==aus8_figures.cross.lon(1,sp)),...
            aus8_currents.lat_v_SBC_north(...
            aus8_currents.lon_u_ALLC==aus8_figures.cross.lon(1,sp))],...
            [aus8_currents.z_mid, aus8_currents.z_top] ...
            , 'g', ...
            'linewidth',line_width)
    end
    
    if aus8_figures.cross.lon(1,sp)~=147
        plot([...
            aus8_currents.lat_v_SBC_south(...
            aus8_currents.lon_u_ALLC==aus8_figures.cross.lon(1,sp)),...
            aus8_currents.lat_v_SBC_north(...
            aus8_currents.lon_u_ALLC==aus8_figures.cross.lon(1,sp))],...
            [aus8_currents.z_mid, aus8_currents.z_mid] ...
            , 'b', ...
            'linewidth',line_width)
    else
        plot([...
            aus8_currents.lat_v_SBC_south(end),...
            aus8_currents.lat_v_SBC_north(end)],...
            [aus8_currents.z_mid, aus8_currents.z_mid] ...
            , 'b', ...
            'linewidth',line_width)
    end
    
    if aus8_figures.cross.lon(1,sp)~=147
        plot([...
            aus8_currents.lat_v_SBC_south(...
            aus8_currents.lon_u_ALLC==aus8_figures.cross.lon(1,sp)),...
            aus8_currents.lat_v_SBC_south(...
            aus8_currents.lon_u_ALLC==aus8_figures.cross.lon(1,sp))],...
            [aus8_currents.z_mid, aus8_currents.z_top] ...
            , 'r', ...
            'linewidth',line_width)
    else
        plot([...
            aus8_currents.lat_v_SBC_south(end),...
            aus8_currents.lat_v_SBC_south(end)],...
            [aus8_currents.z_mid, aus8_currents.z_top] ...
            , 'r', ...
            'linewidth',line_width)
    end
    
    if aus8_figures.cross.lon(1,sp)~=147
        plot([...
            aus8_currents.lat_v_DRC_south(...
            aus8_currents.lon_u_ALLC==aus8_figures.cross.lon(1,sp)),...
            aus8_currents.lat_v_DRC_south(...
            aus8_currents.lon_u_ALLC==aus8_figures.cross.lon(1,sp))],...
            [aus8_currents.z_mid, aus8_currents.z_top] ...
            , ...
            '-m', ...
            'linewidth',line_width)
    else
        plot([...
            aus8_currents.lat_v_DRC_south(end),...
            aus8_currents.lat_v_DRC_south(end)],...
            [aus8_currents.z_mid, aus8_currents.z_top] ...
            , ...
            '-m', ...
            'linewidth',line_width)
    end
    
    if aus8_figures.cross.lon(1,sp)~=147
        plot([...
            aus8_currents.lat_v_DRC_south(...
            aus8_currents.lon_u_ALLC==aus8_figures.cross.lon(1,sp)),...
            aus8_currents.lat_v_DRC_south(...
            aus8_currents.lon_u_ALLC==aus8_figures.cross.lon(1,sp))],...
            [aus8_currents.z_bot, aus8_currents.z_mid] ...
            , ...
            '-m', ...
            'linewidth',line_width)
    else
        plot([...
            aus8_currents.lat_v_DRC_south(end),...
            aus8_currents.lat_v_DRC_south(end)],...
            [aus8_currents.z_bot, aus8_currents.z_mid] ...
            , ...
            '-m', ...
            'linewidth',line_width)
    end
    
    if aus8_figures.cross.lon(1,sp)~=147
        plot([...
            aus8_currents.lat_v_DRC_south(...
            aus8_currents.lon_u_ALLC==aus8_figures.cross.lon(1,sp)),...
            aus8_currents.lat_v_SBC_north(...
            aus8_currents.lon_u_ALLC==aus8_figures.cross.lon(1,sp))],...
            [-2000, -2000] ...
            , ...
            'c', ...
            'linewidth',line_width)
    else
        plot([...
            aus8_currents.lat_v_DRC_south(end),...
            aus8_currents.lat_v_SBC_north(end)],...
            [-2000, -2000] ...
            , ...
            'c', ...
            'linewidth',line_width)
    end
    
    h_tit = title(['(' lett(sp) ') $' ...
        num2str(aus8_figures.cross.lon(1,sp)) ...
        '^{\circ}$E: ' aus8_figures.sect_names{sp}], ...
        'fontsize',font_size, 'horizontalalignment','left');
    h_tit.Position(1) = axis_setup{sp}(1);

    grid
    set(ax,'layer','top','color',nan_color,...
        'fontsize',font_size,'tickdir','out', ...
        'ytick', (-2000:200:0), ...
        'yticklabel', -2000:200:0, ...
        'ticklength',fig_tick_length, ...
        'xtick',axis_setup{sp}(1):axis_setup{sp}(2)-1)
    if col_ind(sp) ~= 1
        set(gca,'yticklabel','')
    else
        ylabel('Depth (m)')
    end
    if row_ind(sp) == 2
        xlabel('Latitude')
    end
    
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
            cbar_x, cbar_y]);
        set(get(cbar,'ylabel'),'String','$u$ (cm/s)', ...
            'fontsize',font_size)
        cbar.Label.Interpreter = 'latex';
        
        pointycbar(cbar)
        
    elseif sp == 4
        ax = axes('visible', 'off');
        colormap(ax, cmap);
        cbar = colorbar;
        set(cbar,'ytick',linspace(0.05,0.95,8), ...
            'YAxisLocation','right', ...
            'YTickLabel',[-0.004:0.001:-0.001,0.001:0.001:0.004],...
            'fontsize',font_size,'ticklength',cbar_tick_length)
        set(cbar,'units','centimeters','position', [...
            (marg_l+x_sp*(cm(sp))+gap_w*(cm(sp)-1)+plot_cbar_gap), ...
            (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
            cbar_x, ...
            cbar_y]);
        set(get(cbar,'ylabel'),'String','$w$ (cm/s)', ...
            'fontsize',font_size)
        cbar.Label.Interpreter = 'latex';
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
print(fig, ...
    ['~/Duran2017/SACS/10319442jbhpxfsdfvwy/' scriptname(1:3) ...
    '_fig' num2str(fig_n) '_'], ...
    '-dpng', '-r300')
% close


%%
play(bird_f_path,[1 (get(bird_f_path, 'SampleRate')*3)]);

