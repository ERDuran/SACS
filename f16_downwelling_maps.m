%% fig 2: meridional cross sections of zonal velocity differences
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/KDau_fulw'])
load([data_path 'SACS_data/aus8_figures'])
load([data_path 'SACS_data/SmSan02'])

topo_binavg = SmSan02.topo_binavg;
topo_binavg(topo_binavg==0) = NaN;

%%
depth = aus8_coor.depth;
lon = aus8_coor.lon_v;
lat = aus8_coor.lat_u;
fulw = KDau_fulw.mean(:,:,31);


%% 5) plot maps of U and V SBC
screen_ratio = 0.75;
fig_n = 1;
rowcols = [1 1];
rowcols_size = [12 4.2]/screen_ratio; % cm
margs = [1.0 0.3 1.9 0.7]/screen_ratio; % cm
gaps = [0.4 0.7]/screen_ratio; % cm

magnif = 100;
cmap1_cont = -[100 50 20 10 5 2.5 0.5 0];
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
cmaps_y_label = cmaps_cont;

% lon_min = 115; lon_max = 147; lat_min = -47; lat_max = -32;
lon_min = 110; lon_max = 152; lat_min = -48; lat_max = -31;

x_chc = {lon};
x_ind = [1];
y_chc = {lat};

w_conv = 10^4;
data = {fulw*magnif*w_conv};

title_chc = ...
    {'{\boldmath{$V_{up}$}} (arrows) and {$U_{up}$} (shadings)'};

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
        
    colormap(ax, cmaps_custom{sp});
    pcolor(x{sp}, y{sp}, data{sp})
    axis(axis_setup{sp})
    shading interp
    caxis([minmax{sp}(1) minmax{sp}(2)]);
    hold on
    
%     contour(aus8_coor.lon,aus8_coor.lat,topo_binavg, [-500 -500], '--')
    contour(aus8_coor.lon,aus8_coor.lat,topo_binavg, [-1000 -1000], 'k-')
    
    for n = 1 : length(aus8_figures.cross.lon(1,:))
        plot(aus8_figures.cross.lon(1:2,n), ...
            aus8_figures.cross.lat(1:2,n), '--k', 'linewidth',1.2)
    end
    
    
    h_tit = title(['(' lett(sp) ') Vertical velocity at z = -250 m depth in MOM01.'], ...
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
        ax = axes('visible', 'off');
        colormap(ax, cmaps);
        cbar = colorbar('horizontal');
        set(cbar,'xtick',cmaps_linspace, ...
            'XTickLabel',...
            [cmaps_y_label(1:7), 0, cmaps_y_label(9:end)],...
            'fontsize',font_size,'ticklength',cbar_tick_length)
        set(cbar,'units','centimeters','position', [...
            marg_l, ...
            marg_b+plot_cbar_gap, ...
            cbar_x, ...
            cbar_y]);
        
        set(get(cbar,'xlabel'), ...
            'String','$w$ ($10^{-4}$ cm/s)', ...
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
    '_fig' num2str(fig_n) '_'], ...
    '-dpng', '-r300')
% print(fig, ...
%     ['~/Duran2017/SACS/10319442jbhpxfsdfvwy/' scriptname(1:3) ...
%     '_fig' num2str(fig_n) '_'], ...
%     '-dpng', '-r300')
% close


%%
play(bird_f_path,[1 (get(bird_f_path, 'SampleRate')*3)]);

