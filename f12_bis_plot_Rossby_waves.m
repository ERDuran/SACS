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

[kds_var, ~, kds_att] = ...
        nc2mat([data_path ...
        'KDS75_cp/sea_level_7yrs.nc'], ...
        'ALL');
    
lon_ori = kds_var.xt_ocean +360;
lat_ori = kds_var.yt_ocean;
time = kds_var.time;
v_orig = kds_var.sea_level;
size(v_orig)


%% 1st interp: Latitude
lat_hr = -42;
v_temp_lat = NaN(length(time), length(lon_ori));
for j = 1 : length(lon_ori)
    for t = 1 : length(time)
        v_temp_lat(t,j) = ...
            interp1(lat_ori, squeeze(v_orig(j,:,t)), lat_hr);
        
    end
    disp(['lon = ' num2str(lon_ori(j))])
end
            

%% 3rd interp: Lon
lon_hr = lon;
V_hr = NaN(length(time), length(lon_hr));
for t = 1 : length(time)
    V_hr(t,:) = interp1(lon_ori, squeeze(v_temp_lat(t,:)), lon_hr);
    
end

V_hr(:,290:320) = NaN;


%% 5) plot maps of U and V SBC
screen_ratio = 0.75;
fig_n = 1;
rowcols = [1 1];
rowcols_size = [13 20]/screen_ratio/2; % cm
margs = [1.4 0.2 1.8 1]/screen_ratio; % cm
gaps = [0.4 0.8]/screen_ratio; % cm
plot_cbar_gap = 1/screen_ratio;
cbar_x = rowcols_size(1);
cbar_y = 0.2/screen_ratio;

magnif = 1000;
cmap1_cont = -[0.3 0.2 0.15 0.1 0.075 0.05 0.025 0.01 0.001 0]*magnif;
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

for sp = 1 : rowcols(1)*rowcols(2)
    minmax{sp} = [cmaps_cont(1) cmaps_cont(end)];
    cmaps_custom{sp} = cmapcust(cmaps,cmaps_cont);

end

Seasons = {'Summer', 'Autumn', 'Winter', 'Spring'};
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
    
    colormap(ax, cmaps_custom{sp});
    pcolor(lon,1:28,V_hr*magnif)
    shading interp
    caxis([minmax{sp}(1) minmax{sp}(2)]);
    axis([110 152 1 28])
    
    %title_chc = ...
    %{'{\boldmath{$V_{up}$}} (arrows) and {$U_{up}$} (shadings)', ...
    %'{\boldmath{$V_{low}$}} (arrows) and {$U_{low}$} (shadings)'};
    
    h_tit = title({'MOM01-75z seasonality of sea level anomaly along $40^{\circ}S$'}, ...
    'horizontalalignment','left', 'fontsize',font_size);
    h_tit.Position(1) = 110;
    grid
    set(ax,'layer','top','color',nan_color,...
        'fontsize',font_size,'tickdir','out', ...
        'ticklength',fig_tick_length, ...
        'xtick', 110:4:152, 'ytick', 1:28, 'yticklabel', Seasons)
    if row_ind(sp) ~= rowN
        set(gca,'xticklabel','')
    else
        xlabel('Longitude')
    end
    if col_ind(sp) ~= 1
        set(gca,'yticklabel','')
    else
        ylabel('Season')
    end
    
    hold on
    plot([110 152], [2 2]+2, 'k--')
    plot([110 152], [6 6]+2, 'k--')
    plot([110 152], [10 10]+2, 'k--')
    plot([110 152], [14 14]+2, 'k--')
    plot([110 152], [18 18]+2, 'k--')
    plot([110 152], [22 22]+2, 'k--')
    plot([110 152], [26 26]+2, 'k--')
    
    gapi = 1.8;
    if sp == rowN*colN
        ax = axes('visible', 'off');
        colormap(ax, cmaps);
        cbar = colorbar('horizontal');
        set(cbar,'ytick',cmaps_linspace, ...
            'YAxisLocation','right','YTickLabel',cmaps_y_label*magnif,...
            'fontsize',font_size,'ticklength',cbar_tick_length)
        set(cbar,'units','centimeters','position', [...
            (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1))-gapi, ...
            (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)-plot_cbar_gap), ...
            cbar_x+gapi, ...
            cbar_y]);
        set(get(cbar,'xlabel'),'String','$\eta$ ($mm$)', ...
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
    [figures_path mfilename '/f12_plot_Rossby_waves/f12' ...
    '_fig' num2str(fig_n) '_'], ...
    '-dpng', '-r300')
print(fig, ...
    ['~/Duran2017/SACS/10319442jbhpxfsdfvwy/f12' ...
    '_fig' num2str(fig_n) '_'], ...
    '-dpng', '-r300')
close


%%
play(bird_f_path,[1 (get(bird_f_path, 'SampleRate')*3)]);

