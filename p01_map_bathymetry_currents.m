%% fig 2: map bathymetry and landmarks. In black and white
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/SmSan02'])

lat = aus8_coor.lat;
lon = aus8_coor.lon;
topo_binavg = SmSan02.topo_binavg;
topo_binavg(topo_binavg==0) = NaN;


%% plot map
screen_ratio = 0.75;
fig_n = 1;
rowcols = [1 1];
rowcols_size = [13.8 6]/screen_ratio; % cm
margs = [1.4 2 1.2 0.2]/screen_ratio; % cm
gaps = [1 1]/screen_ratio; % cm
plot_cbar_gap = 0.2/screen_ratio;
cbar_x = 0.2/screen_ratio;
cbar_y = rowcols_size(2);

for sp = 1 : rowcols(1)*rowcols(2)
    axis_setup{sp} = ...
        [114 148 ...
        -46 -31];
    x{sp} = aus8_coor.lon;
    y{sp} = aus8_coor.lat;
end
    
data = {topo_binavg};
minmax = {[-2000 0]};

cmap1_cont = ...
    [-2000 : 200 : -200];
cmap2_cont = ...
    [-200 : 10 : 0];
lvl_cmap1 = length(cmap1_cont)-1;
lvl_cmap2 = length(cmap2_cont)-1;
cmap1 = othercolor('BuPu9', lvl_cmap1);
cmap2 = flipud(othercolor('Paired8', lvl_cmap2));
cmaps = [cmap1; cmap2];
cmaps_cont = [cmap1_cont cmap2_cont(2:end)];
cmaps_cont_length = length(cmaps_cont);
cmaps_linspace = linspace(0,1,cmaps_cont_length);
cmaps_custom = {cmapcust(cmaps,cmaps_cont)};

font_size = 12*screen_ratio;
nan_color = [0 0 0];
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
    'position',[0 0 fig_x fig_y], ...
    'paperposition',[0 0 fig_x fig_y]*screen_ratio);

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
    
    % sections
    merid_Xsect_lat = [ ...
        -34; %1
        -35; %2
        -35; %3
        -35; %4
        -34; %5
        -34; %6
        -34; %7
        -34; %8
        -34; %9
        -34; %10
        -33; %11
        -33; %12
        -32; %13
        -32; %14
        -32; %15
        -32; %16
        -32; %17
        -32; %18
        -33; %19
        -33; %20
        -34; %21
        -35; %22
        -36; %23
        -36; %24
        -37; %25
        -37; %26
        -38; %27
        -38; %28
        -39; %29
        -40; %30
        -42; %31
        -43; %32
        -43; %33
        ]';
    merid_Xsect_lat(2,:) = merid_Xsect_lat(1,:) - 5;
    merid_Xsect_lon = [115:147; 115:147];
    
    aus8_figures.sect_names = {...
        'C. Leeuwin', 'C. Pasley', 'west GAB', 'east GAB', ...
        'C. Carnot', 'Portland', 'King Is.', 'South East C.'};
        
    lon_n = [1 9 14 18 22 27 30 33];
    
    small_arr = 3;
    big_ar_head = 12;
    FC_south = 0.65;
    % FC
    FC_color = [1 0 0];
    contour(x{sp}(62:268), y{sp}-FC_south, data{sp}(:,62:268), ...
        [-700 -700], ...
        'color',FC_color, 'linewidth', 1.6)
    contour(x{sp}(272:313)-0.6, y{sp}-FC_south+0.1, ...
        data{sp}(:,272:313), ...
        [-700 -700], ...
        'color',FC_color, 'linewidth', 1.6)
    contour(x{sp}(308:313), y{sp}-FC_south+0.05, data{sp}(:,308:313), ...
        [-700 -700], ...
        'color',FC_color, 'linewidth', 1.6)
    arrow([115.1,-35.2-FC_south], [115,-35.2-FC_south], ...
        big_ar_head, ...
        'Facecolor', FC_color, 'edgecolor', FC_color)
    text(126.5,-34.5-FC_south, 'Flinders Current', 'fontsize', font_size, ...
        'color', FC_color)
    arrow([146,-44.0-FC_south], [143,-43.8-FC_south], ...
        big_ar_head, ...
        'Facecolor', FC_color, 'edgecolor', FC_color, ...
        'linewidth', 1.6)
    text(141,-44.5-FC_south, 'Tasman Leakage', 'fontsize', font_size, ...
        'color', FC_color)
    
    % LCE
    SBC_color = [0 0 1];
    LC_south = 0.2;
    contour(x{sp}(57:127), y{sp}-LC_south, data{sp}(:,57:127), ...
        [-700 -700], ...
        'color',SBC_color, 'linewidth', 1.3)
    arrow([124,-34.49-LC_south], [124.1,-34.44-LC_south], ...
        big_ar_head-2, ...
        'Facecolor', SBC_color, 'edgecolor', SBC_color)
    text(121,-37.6-LC_south, {'Leeuwin Current','Extension'}, ...
        'fontsize', font_size, ...
        'color', SBC_color,'horizontalalignment','right')
    arrow([121,-37.6-LC_south], [122,-34.8], small_arr, ...
        'Facecolor', [0 0 0], 'edgecolor', [0 0 0])

    % SAC
    SAC_south = 0.275;
    contour(x{sp}(132:263), y{sp}-SAC_south, data{sp}(:,132:263), ...
        [-700 -700], ...
        'color',SBC_color, 'linewidth', 1.3)
    arrow([141.0,-38.6-SAC_south], [141.1,-38.65-SAC_south], ...
        big_ar_head-2, ...
        'Facecolor', SBC_color, 'edgecolor', SBC_color)
    text(134,-36.5-SAC_south, 'South Australian Current', 'fontsize', font_size, ...
        'color', SBC_color,'horizontalalignment','right')
    arrow([134,-36.5-SAC_south], [135.25,-36], small_arr, ...
        'Facecolor', [0 0 0], 'edgecolor', [0 0 0])
    
    % ZC
    ZC_south = 0.35;
    contour(x{sp}(266:310), y{sp}-ZC_south, data{sp}(:,266:310), ...
        [-700 -700], ...
        'color',SBC_color, 'linewidth', 1.3)
    arrow([146.9,-44.21-ZC_south], [147,-44.24-ZC_south], ...
        big_ar_head-2, ...
        'Facecolor', SBC_color, 'edgecolor', SBC_color)
    text(143,-42-ZC_south, 'Zeehan Current', 'fontsize', font_size, ...
        'color', SBC_color,'horizontalalignment','right')
    arrow([143,-42-ZC_south], [144.5,-42-ZC_south], small_arr, ...
        'Facecolor', [0 0 0], 'edgecolor', [0 0 0])

    grid
    set(ax,'layer','top','color',nan_color,...
        'fontsize',font_size,'tickdir','out','ticklength',fig_tick_length)
    if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), 
    else, set(gca,'xtick',108:2:152,'ytick',-50:2:-30), end
    if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end
end


xlabel('Longitude', 'fontsize', font_size)
ylabel('Latitude', 'fontsize', font_size)

ax = axes('visible', 'off');
colormap(ax, cmaps);
cbar = colorbar;
set(cbar,'ytick',cmaps_linspace(2:2:end), ...
    'YTickLabel',cmaps_cont(2:2:end),'fontsize',font_size,...
    'ticklength',cbar_tick_length)
set(cbar,'units','centimeters', ...
    'position', [...
    (marg_l+x_sp*(cm(sp))+gap_w*(cm(sp)-1)+plot_cbar_gap), ...
    (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
    cbar_x, ...
    cbar_y]);

set(get(cbar,'ylabel'),'String','Depth (m)','fontsize',font_size)
cbar.Label.Interpreter = 'latex';

outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
print(fig, ...
    [figures_path mfilename '/' scriptname(1:3) ...
    '_fig' num2str(fig_n) '_'], ...
    '-dpng', '-r300')
% close

%%
save([data_path 'SACS_data/aus8_figures'], 'aus8_figures')
disp('aus8_figures DONE')
play(bird_f_path,[1 (get(bird_f_path, 'SampleRate')*3)]);

