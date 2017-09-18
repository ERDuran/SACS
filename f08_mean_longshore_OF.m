%%
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_currents'])
load([data_path 'SACS_data/aus8_figures'])
load([data_path 'SACS_data/aus8_u_g_prime'])
load([data_path 'SACS_data/aus8_v_g_prime'])
load([data_path 'SACS_data/KDau_fulu'])
load([data_path 'SACS_data/KDau_fulv'])

a = aus8_coor.a;
pi180 = aus8_coor.pi180;
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;
lat_v_SBC_north = aus8_currents.lat_v_SBC_north;
lat_v_SBC_south = aus8_currents.lat_v_SBC_south;
lat_v_DRC_north = aus8_currents.lat_v_DRC_north;
lat_v_DRC_south = aus8_currents.lat_v_DRC_south;
lon_u_ALLC = aus8_currents.lon_u_ALLC;
Months = aus8_coor.Months;
Months{13} = 'mean';
depth_mid = aus8_coor.depth_mid;

u_g_prime = aus8_u_g_prime.mean;
v_g_prime = aus8_v_g_prime.mean;


%% DRC Vts up
vt_lon_v = lon_v(57:312);
vt_g_prime_DRC = zeros(length(depth_mid), length(vt_lon_v));

% in between cases
for jj = 1 : length(vt_lon_v)
    lon_v_ind = find(lon_v == vt_lon_v(jj));
    
    %
    lat_ind_V_Vts_prime_up_DRC = lat_v == lat_v_DRC_south(jj);
    V_Vts_prime_up_DRC_now = ...
        v_g_prime(lat_ind_V_Vts_prime_up_DRC, lon_v_ind, :);
    V_Vts_prime_up_DRC = V_Vts_prime_up_DRC_now;
    
    vt_g_prime_DRC(:,jj) = squeeze(V_Vts_prime_up_DRC);
end


%% 5) plot maps of U and V SBC
depth_mid_for_plot = aus8_coor.depth_mid;
depth_mid_for_plot(end) = -2000;

screen_ratio = 0.75;
fig_n = 1;
rowcols = [1 1];
rowcols_size = [14 6]/screen_ratio/2; % cm
margs = [0.6 0.2 1.2 0.6]/screen_ratio; % cm
gaps = [0.4 0.8]/screen_ratio; % cm
plot_cbar_gap = 0.7/screen_ratio;
cbar_x = rowcols_size(1);
cbar_y = 0.2/screen_ratio;

magnif = 1000;
cmap1_cont = -[20 10 8 6 4 2 1 0.5 0.1 0]*magnif;
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


minmax = [cmaps_cont(1) cmaps_cont(end)];
cmaps_custom = cmapcust(cmaps,cmaps_cont);

x = vt_lon_v;
y = depth_mid_for_plot;

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
    
    colormap(ax, cmaps_custom);
    pcolor(x, y, vt_g_prime_DRC*100*magnif)
    axis([115 147 -2000 0])
    shading interp
    caxis([minmax(1) minmax(2)]);
    
    hold on
    plot([115 147], [-250 -250], '-b', 'linewidth', 0.3)
    
    
    h_tit = title(['(' lett(sp) ') ' ...
        ' Mean ZD CARS $V_g''$ ($cm/s$) along $y_{OF}$'], ...
    'horizontalalignment','left', 'fontsize',font_size);
    h_tit.Position(1) = 115;
    grid
    set(ax,'layer','top','color',nan_color,...
        'fontsize',font_size,'tickdir','out', ...
        'ticklength',fig_tick_length, ...
        'xtick', 115:2:147, ...
        'ytick', -2000:200:0)
%     if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
%     if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end
    
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

