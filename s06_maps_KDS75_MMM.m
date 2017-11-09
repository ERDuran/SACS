%%
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_currents'])
load([data_path 'SACS_data/KDau_fcrt'])
load([data_path 'SACS_data/aus8_figures'])

MTH = aus8_coor.MTH;
lat = aus8_coor.lat;
lon = aus8_coor.lon;
depth_mid = aus8_coor.depth_mid;
depth = aus8_coor.depth;
depth_thkn = aus8_coor.depth_thkn;
SEASONS = {'Summer', 'Autumn', 'Winter', 'Spring'};


%%
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;

for t = 1 : 4
    fulu_ztop_to_zmid.(MTH{t}) = ...
        KDau_fcrt.MMM.ztop_to_zmid.fulu.(MTH{t});
    fulv_ztop_to_zmid.(MTH{t}) = ...
        KDau_fcrt.MMM.ztop_to_zmid.fulv.(MTH{t});
    fulv_ztop_to_zmid_interp2.(MTH{t}) = interp2(...
        lon_v, lat_v, fulv_ztop_to_zmid.(MTH{t}), lon_u, lat_u);
    
    fulu_zmid_to_zbot.(MTH{t}) = ...
        KDau_fcrt.MMM.zmid_to_zbot.fulu.(MTH{t});
    fulv_zmid_to_zbot.(MTH{t}) = ...
        KDau_fcrt.MMM.zmid_to_zbot.fulv.(MTH{t});
    fulv_zmid_to_zbot_interp2.(MTH{t}) = interp2(...
        lon_v, lat_v, fulv_zmid_to_zbot.(MTH{t}), lon_u, lat_u);
end

z_top = aus8_currents.z_top;
z_mid = aus8_currents.z_mid;
z_bot = aus8_currents.z_bot;


%% 5) plot maps of U and V SBC
screen_ratio = 0.75;
fig_n = 1;
rowcols = [2 1];
rowcols_size = [14 6]/screen_ratio/2; % cm
margs = [0.6 0.2 1.2 0.6]/screen_ratio; % cm
gaps = [0.4 0.8]/screen_ratio; % cm
plot_cbar_gap = 0.7/screen_ratio;
cbar_x = rowcols_size(1);
cbar_y = 0.2/screen_ratio;

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
lon_min = 108; lon_max = 153; lat_min = -50; lat_max = -30;

x_chc = {aus8_coor.lon_u, aus8_coor.lon_v};
x_ind = [1 1];
y_chc = {aus8_coor.lat_u, aus8_coor.lat_v};

for t = 1 : 4
    data{t} = fulu_ztop_to_zmid.(MTH{t})*magnif;
    v_data{t} = fulv_ztop_to_zmid_interp2.(MTH{t})*magnif;
end

title_chc = {'U', 'U', 'U', 'U'};
z1_chc = {z_top, z_mid};
z2_chc = {z_mid, z_bot};

for sp = 1 : rowcols(1)*rowcols(2)
    minmax{sp} = [cmaps_cont(1) cmaps_cont(end)];
    cmaps_custom{sp} = cmapcust(cmaps,cmaps_cont);
    
    axis_setup{sp} = [lon_min lon_max lat_min lat_max];
    x{sp} = x_chc{x_ind(1)};
    y{sp} = y_chc{x_ind(1)};

end

lett = 'a':'z';
font_size = 10*screen_ratio;
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
cbar_tick_length = norm_length;

set(fig,'units','centimeters','paperunits','centimeters', ...
    'inverthardcopy','off','color',fig_color,...
    'paperposition',[0 0 fig_x fig_y]*screen_ratio,...
    'position',[0 0 fig_x fig_y]);

for sp = 1 : 2
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
    
    for n = 1 : length(aus8_figures.cross.lon(1,:))
        plot(aus8_figures.cross.lon(1:2,n), ...
            aus8_figures.cross.lat(1:2,n), '--k')
    end
    
    
    h_tit = title(SEASONS{sp}, ...
    'horizontalalignment','left', 'fontsize',font_size);
    h_tit.Position(1) = axis_setup{sp}(1);
    grid
    set(ax,'layer','top','color',nan_color,...
        'fontsize',font_size,'tickdir','out', ...
        'ticklength',fig_tick_length, ...
        'xtick', lon_min:4:lon_max, 'ytick', lat_min:4:lat_max)
    if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
    if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end
    
    if sp == rowN*colN
        ax = axes('visible', 'off');
        colormap(ax, cmaps);
        cbar = colorbar('horizontal');
        set(cbar,'ytick',cmaps_linspace, ...
            'YAxisLocation','right','YTickLabel',cmaps_y_label,...
            'fontsize',font_size,'ticklength',cbar_tick_length)
        set(cbar,'units','centimeters','position', [...
            (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)), ...
            (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)-plot_cbar_gap), ...
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

