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
fig_n = 1;
rowcols = [1 1];
rowcols_size = [25 15]; % cm
margs = [1.2 2.5 1.2 1.2]; % cm
gaps = [1 1]; % cm
plot_cbar_gap = 0.5;
cbar_x = 0.3;
cbar_y = rowcols_size(2);

for sp = 1 : rowcols(1)*rowcols(2)
    axis_setup{sp} = ...
        [aus8_coor.lon(1) aus8_coor.lon(end) ...
        aus8_coor.lat(end) aus8_coor.lat(1)];
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

font_size = 13;
fig_color = [0 0 0];


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
set(fig,'units','centimeters',...
    'position',[0 0 ...
    (marg_l+colN*x_sp+gap_w*(colN-1)+marg_r) ...
    (marg_b+rowN*y_sp+gap_h*(rowN-1)+marg_t)]);

for sp = 1 : rowN*colN
    ax = axes('Units','centimeters', ...
            'Position',[...
            (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)), ...
            (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
            x_sp, ...
            y_sp]);
        
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
    
    for n = 1 : length(lon_n)
        aus8_figures.cross.lon(1:2,n) = merid_Xsect_lon(:,lon_n(n));
        aus8_figures.cross.lat(1:2,n) = merid_Xsect_lat(:,lon_n(n));
        plot(merid_Xsect_lon(:,lon_n(n)), ...
            merid_Xsect_lat(:,lon_n(n)), '--k')
    end
    
    big_ar_head = 18;
    FC_south = 1.2;
    % FC
    FC_color = [1 0 0];
    contour(x{sp}(66:268), y{sp}-FC_south, data{sp}(:,66:268), ...
        [-700 -700], ...
        'color',FC_color, 'linewidth', 3)
    contour(x{sp}(272:312)-0.5, y{sp}-FC_south+0.1, data{sp}(:,272:312), ...
        [-700 -700], ...
        'color',FC_color, 'linewidth', 3)
    contour(x{sp}(308:313), y{sp}-FC_south+0.05, data{sp}(:,308:313), ...
        [-700 -700], ...
        'color',FC_color, 'linewidth', 3)
    arrow([115.1,-35.4-FC_south], [115,-35.4-FC_south], ...
        big_ar_head, ...
        'Facecolor', FC_color, 'edgecolor', FC_color)
    text(133,-35.7-FC_south, 'FC', 'fontsize', font_size, ...
        'color', FC_color)
    
    % LC
    SBC_color = [0 0 1];
    LC_south = 0.2;
    contour(x{sp}(57:175), y{sp}-LC_south, data{sp}(:,57:175), ...
        [-700 -700], ...
        'color',SBC_color, 'linewidth', 2)
    arrow([130,-33.6-LC_south], [130.1,-33.6-LC_south], ...
        big_ar_head-4, ...
        'Facecolor', SBC_color, 'edgecolor', SBC_color)
    text(121,-35.1-LC_south, 'LC', 'fontsize', font_size, ...
        'color', SBC_color)
    
    % SAC
    SAC_south = 0.275;
    contour(x{sp}(183:272), y{sp}-SAC_south, data{sp}(:,183:272), ...
        [-700 -700], ...
        'color',SBC_color, 'linewidth', 2)
    arrow([142.1,-39.15-SAC_south], [142.2,-39.2-SAC_south], ...
        big_ar_head-4, ...
        'Facecolor', SBC_color, 'edgecolor', SBC_color)
    text(137,-37.2-SAC_south, 'SAC', 'fontsize', font_size, ...
        'color', SBC_color, 'rotation', -35)
    
    % ZC
    ZC_south = 0.5;
    contour(x{sp}(278:306), y{sp}-ZC_south, data{sp}(:,278:306), ...
        [-700 -700], ...
        'color',SBC_color, 'linewidth', 2)
    arrow([146.9,-44.21-ZC_south], [147,-44.24-ZC_south], ...
        big_ar_head-4, ...
        'Facecolor', SBC_color, 'edgecolor', SBC_color)
    text(143.1,-41-ZC_south, 'ZC', 'fontsize', font_size, ...
        'color', SBC_color, ...
        'rotation', -55)
    
    % CC
    CC1 = [127, -32.7];
    arrow([CC1(1), CC1(2)], [CC1(1)+4, CC1(2)], ...
        big_ar_head-8, 'linewidth', 2, ...
        'Facecolor', [1 1 0], 'edgecolor', [1 1 0])
    arrow([CC1(1)+4, CC1(2)], [CC1(1), CC1(2)], ...
        big_ar_head-8, 'linewidth', 2, ...
        'Facecolor', [1 1 0], 'edgecolor', [1 1 0])
    CC2 = [132.7, -32.8];
    arrow([CC2(1), CC2(2)], [CC2(1)+1.8, CC2(2)-1.8], ...
        big_ar_head-8, 'linewidth', 2, ...
        'Facecolor', [1 1 0], 'edgecolor', [1 1 0])
    arrow([CC2(1)+1.8, CC2(2)-1.8], [CC2(1), CC2(2)], ...
        big_ar_head-8, 'linewidth', 2, ...
        'Facecolor', [1 1 0], 'edgecolor', [1 1 0])
    CC3 = [149, -31]; % reference arrow
    arrow([CC3(1), CC3(2)], [CC3(1)+2.3, CC3(2)], ...
        big_ar_head-8, 'linewidth', 2, ...
        'Facecolor', [1 1 0], 'edgecolor', [1 1 0])
    arrow([CC3(1)+2.3, CC3(2)], [CC3(1), CC3(2)], ...
        big_ar_head-8, 'linewidth', 2, ...
        'Facecolor', [1 1 0], 'edgecolor', [1 1 0])
    text(CC3(1)+0.5, CC3(2)-0.7, 'CC', 'fontsize', font_size, ...
        'color', [1 1 0])
    CC4 = [145, -39.9];
    arrow([CC4(1), CC4(2)], [CC4(1)+2.3, CC4(2)], ...
        big_ar_head-8, 'linewidth', 2, ...
        'Facecolor', [1 1 0], 'edgecolor', [1 1 0])
    arrow([CC4(1)+2.3, CC4(2)], [CC4(1), CC4(2)], ...
        big_ar_head-8, 'linewidth', 2, ...
        'Facecolor', [1 1 0], 'edgecolor', [1 1 0])

    % OF
    OF1 = [120, -36.5];
    arrow([OF1(1), OF1(2)-2], [OF1(1), OF1(2)], ...
        big_ar_head-8, 'linewidth', 2, ...
        'Facecolor', [1 0.5 0], 'edgecolor', [1 0.5 0])
    text(OF1(1)+0.5, OF1(2)-1, 'OF', 'fontsize', font_size, ...
        'color', [1 0.5 0])
    OF2 = [138, -39];
    arrow([OF2(1), OF2(2)-2], [OF2(1), OF2(2)], ...
        big_ar_head-8, 'linewidth', 2, ...
        'Facecolor', [1 0.5 0], 'edgecolor', [1 0.5 0])
    OF3 = [130, -35.5];
    arrow([OF3(1), OF3(2)-2], [OF3(1), OF3(2)], ...
        big_ar_head-8, 'linewidth', 2, ...
        'Facecolor', [1 0.5 0], 'edgecolor', [1 0.5 0])

    grid
    set(gca,'layer','top','color',fig_color,...
        'fontsize',font_size,'tickdir','out')
    if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), 
    else set(gca,'xtick',108:2:152), end
    if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end
end


% 8) place landmarks
small_arr = 6;
% Cape Leeuwin
arrow([117,-32], [115,-34.15], small_arr, ...
    'Facecolor', [1 1 1], 'edgecolor', [1 1 1])
text(117,-32,'Cape Leeuwin', 'fontsize', font_size, ...
    'color', [1 1 1])
% Albany
arrow([119,-33], [118,-34.9], small_arr, ...
    'Facecolor', [1 1 1], 'edgecolor', [1 1 1])
text(119,-33,'Albany', 'fontsize', font_size, ...
    'color', [1 1 1])
% Cape Pasley
arrow([125,-31], [123.4,-33.85], small_arr, ...
    'Facecolor', [1 1 1], 'edgecolor', [1 1 1])
text(125,-31,'Cape Pasley', 'fontsize', font_size, ...
    'color', [1 1 1])
% Great Australian Bight
arrow([131,-30.5], [130,-32.3], small_arr, ...
    'Facecolor', [1 1 1], 'edgecolor', [1 1 1])
text(131,-30.5,'Great Australian Bight', 'fontsize', font_size, ...
    'color', [1 1 1])
% Cape Carnot
arrow([137,-32], [135.7,-34.7], small_arr, ...
    'Facecolor', [1 1 1], 'edgecolor', [1 1 1])
text(137,-32,sprintf('Cape\nCarnot'), 'fontsize', font_size, ...
    'color', [1 1 1])
% Spencer Gulf
arrow([142,-31.5], [136.9,-34.2], small_arr, ...
    'Facecolor', [1 1 1], 'edgecolor', [1 1 1])
text(142,-31.5,'Spencer Gulf', 'fontsize', font_size, ...
    'color', [1 1 1])
% St. Vincent Gulf
arrow([142,-32.5], [138.1,-35], small_arr, ...
    'Facecolor', [1 1 1], 'edgecolor', [1 1 1])
text(142,-32.5,'St. Vincent Gulf', 'fontsize', font_size, ...
    'color', [1 1 1])
% Kangaroo Is.
arrow([142,-34], [137,-35.85], small_arr, ...
    'Facecolor', [1 1 1], 'edgecolor', [1 1 1])
text(142,-34,sprintf('Kangaroo\nIs.'), 'fontsize', font_size, ...
    'color', [1 1 1])
% Cape Jaffa
arrow([140.5,-36], [139.7,-37], small_arr, ...
    'Facecolor', [1 1 1], 'edgecolor', [1 1 1])
text(140.5,-36,sprintf('Cape\nJaffa'), 'fontsize', font_size, ...
    'color', [1 1 1])
% Portland
arrow([145,-35], [141.3,-38], small_arr, ...
    'Facecolor', [1 1 1], 'edgecolor', [1 1 1])
text(145,-35,'Portland', 'fontsize', font_size, ...
    'color', [1 1 1])
% Cape Otway
arrow([145,-36], [143.5,-38.75], small_arr, ...
    'Facecolor', [1 1 1], 'edgecolor', [1 1 1])
text(145,-36,'Cape Otway', 'fontsize', font_size, ...
    'color', [1 1 1])
% Bass Strait
arrow([147,-37.25], [145.6,-39.4], small_arr, ...
    'Facecolor', [1 1 1], 'edgecolor', [1 1 1])
text(147,-37.25,sprintf('Bass\nStrait'), 'fontsize', font_size, ...
    'color', [1 1 1])
% King Is.
arrow([149.3,-39], [144,-39.8], small_arr, ...
    'Facecolor', [0 0 0], 'edgecolor', [0 0 0])
text(149.3,-39,sprintf('King Is.'), 'fontsize', font_size, ...
    'color', [0 0 0])
% Schouten Is.
arrow([149.3,-40.8], [148.2 -42.3], small_arr, ...
    'Facecolor', [0 0 0], 'edgecolor', [0 0 0])
text(149.3,-40.8,sprintf('Schouten\nIs.'), 'fontsize', font_size, ...
    'color', [0 0 0])
% South East Cape
arrow([148,-45.4], [146.7,-43.62], small_arr, ...
    'Facecolor', [0 0 0], 'edgecolor', [0 0 0])
text(148,-45.4,sprintf('South East\nCape'), 'fontsize', font_size, ...
    'color', [0 0 0])
% Tasman Rise
text(149.5,-48,sprintf('Tasman\nRise'), 'fontsize', font_size, ...
    'color', [0 0 0])
% Tasman Sea
arrow([152,-43], [153,-43], small_arr, ...
    'Facecolor', [0 0 0], 'edgecolor', [0 0 0])
text(152,-43,sprintf('Tasman\nSea'), 'fontsize', font_size, ...
    'HorizontalAlignment','right','color', [0 0 0])
% South Australian Basin
text(125,-38.5,sprintf('South Australian Basin'), 'fontsize', font_size, ...
    'color', [0 0 0])

ax = axes;
set(gca,'Visible','off')
colormap(ax, cmaps);
cbar = colorbar;
set(cbar,'ytick',cmaps_linspace(2:2:end), ...
    'YTickLabel',cmaps_cont(2:2:end),'fontsize',font_size)
set(cbar,'units','centimeters','position', [...
    (marg_l+x_sp*(cm(sp))+gap_w*(cm(sp)-1)+plot_cbar_gap), ...
    (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
    cbar_x, ...
    cbar_y]);

outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
export_fig(fig, ...
    [figures_path mfilename '/' scriptname(1:3) ...
    '_fig' num2str(fig_n) '_'], ...
    '-transparent', '-m3')
close

%%
save([data_path 'SACS_data/aus8_figures'], 'aus8_figures')
disp('aus8_figures DONE')
play(bird_f_path,[1 (get(bird_f_path, 'SampleRate')*3)]);

