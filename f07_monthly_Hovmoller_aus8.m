%% fig 2: meridional cross sections of zonal velocity differences
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
Months = aus8_coor.Months;

SBC_Ut = NaN(length(Months), length(lon_u));
[SBC_Vtnc, SBC_Vtsc, SBC_Wtc_real] = ...
    deal(NaN(length(Months), length(lon_v)));
DRC_Ut = NaN(length(Months), length(lon_u));
[DRC_Vtnc, DRC_Vtsc, DRC_Wtc_real] = ...
    deal(NaN(length(Months), length(lon_v)));
for t = 1 : 12
    SBC_Ut(t,:) = aus8_currents.SBC_Ut.(Months{t});
    SBC_Vtnc(t,:) = aus8_currents.SBC_Vtnc.(Months{t});
    SBC_Vtsc(t,:) = aus8_currents.SBC_Vtsc.(Months{t});
    SBC_Wtc_real(t,:) = aus8_currents.SBC_Wtc_real.(Months{t});
    DRC_Ut(t,:) = aus8_currents.DRC_Ut.(Months{t});
    DRC_Vtnc(t,:) = aus8_currents.DRC_Vtnc.(Months{t});
    DRC_Vtsc(t,:) = aus8_currents.DRC_Vtsc.(Months{t});
    DRC_Wtc_real(t,:) = aus8_currents.DRC_Wtc_real.(Months{t});    
end


%%
screen_ratio = 0.75;
fig_n = 1;
rowcols = [2 4];
rowcols_size = [3.7 3.7]/screen_ratio; % cm
margs = [0.6 0.2 1.2 0.4]/screen_ratio; % cm
gaps = [0.2 1.9]/screen_ratio; % cm
plot_cbar_gap = 0.8/screen_ratio;
cbar_x = rowcols_size(1);
cbar_y = 0.2/screen_ratio;

magnif = 100;
cmap1_cont = -[4 3 2 1.5 1 0.5 0.25 0.1 0]*magnif;
cmap2_cont = -fliplr(cmap1_cont);
lvl_cmap1 = length(cmap1_cont)-1;
lvl_cmap2 = length(cmap2_cont)-1;
cmap1 = flipud(othercolor('Blues8', lvl_cmap1));
cmap2 = othercolor('Reds8', lvl_cmap2);
cmap1(end,:) = [1 1 1];
cmap2(1,:) = [1 1 1];
cmaps_top = [cmap1; cmap2];
cmaps_cont_top = [cmap1_cont cmap2_cont(2:end)];
cmaps_cont_length = length(cmaps_cont_top);
cmaps_linspace_top = linspace(0,1,cmaps_cont_length);
cmaps_y_label_top = cmaps_cont_top/magnif;

cmap1_cont = -[17.5 15 12.5 10 7.5 5 2.5 1 0.5 0]*magnif;
cmap2_cont = -fliplr(cmap1_cont);
lvl_cmap1 = length(cmap1_cont)-1;
lvl_cmap2 = length(cmap2_cont)-1;
cmap1 = flipud(othercolor('Blues8', lvl_cmap1));
cmap2 = othercolor('Reds8', lvl_cmap2);
cmap1(end,:) = [1 1 1];
cmap2(1,:) = [1 1 1];
cmaps_bot = [cmap1; cmap2];
cmaps_cont_bot = [cmap1_cont cmap2_cont(2:end)];
cmaps_cont_length = length(cmaps_cont_bot);
cmaps_linspace_bot = linspace(0,1,cmaps_cont_length);
cmaps_y_label_bot = cmaps_cont_bot/magnif;

cmaps_tb = {cmaps_top, cmaps_bot};
cmaps_cont_tb = {cmaps_cont_top, cmaps_cont_bot};
cmaps_chc = [1 1 1 1 2 2 2 2];

lon_min = 115; lon_max = 147; lat_min = 1; lat_max = 12;

x_chc = {aus8_coor.lon_u, aus8_coor.lon_v};
x_ind = [1 2 2 2 1 2 2 2 2];

data = {SBC_Ut*magnif, ...
    SBC_Vtsc*magnif, SBC_Wtc_real*magnif, SBC_Vtnc*magnif, ...
    DRC_Ut*magnif, ...
    DRC_Vtnc*magnif, DRC_Wtc_real*magnif, DRC_Vtsc*magnif};

data_label = {'$\mathcal{U}_{SBC}$', ...
    '$\mathcal{V}_{SBC}$', ...
    '$\mathcal{W}_{SBC}$', '$\mathcal{V}_{CC}$', ...
    '$\mathcal{U}_{FC}$', ...
    '$\mathcal{V}_{FC}$', ...
    '$\mathcal{W}_{FC}$', '$\mathcal{V}_{OF}$'};

title_chc = {'U''', 'V''', 'U_{g}''', 'V_{g}'''};

for sp = 1 : rowcols(1)*rowcols(2)
    cmaps{sp} = cmaps_tb{cmaps_chc(sp)};
    cmaps_cont{sp} = cmaps_cont_tb{cmaps_chc(sp)};
    cmaps_custom{sp} = cmapcust(cmaps{sp},cmaps_cont{sp});
    
    minmax{sp} = [cmaps_cont{sp}(1) cmaps_cont{sp}(end)];
    
    axis_setup{sp} = [lon_min lon_max lat_min lat_max];
    x{sp} = x_chc{x_ind(sp)};
    y{sp} = 1:12;

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
if cbar_y > cbar_x*4+gap_w*3, long_side = cbar_y;
else, long_side = cbar_x*4+gap_w*3; end
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
    
    colormap(ax, cmaps_custom{sp});
    pcolor(x{sp}, y{sp}, data{sp})
    axis(axis_setup{sp})
    shading interp
    caxis([minmax{sp}(1) minmax{sp}(2)]);
    hold on
    
    h_tit = title(['(' lett(sp) ') ' data_label{sp}], ...
        'horizontalalignment','left', 'fontsize',font_size);
    h_tit.Position(1) = axis_setup{sp}(1);
    grid
    x_tick = lon_min:2:lon_max;
    x_tick_label = num2cell(x_tick);
    for c = 1:2:17
        x_tick_label{c} = '';
    end
    
    set(ax,'layer','top','color',nan_color,...
        'ticklength',fig_tick_length, 'xticklabelrotation',90, ...
        'fontsize',font_size,'tickdir','out', ...
        'xtick', x_tick, 'xticklabel', x_tick_label, ...
        'ytick', lat_min:1:lat_max, 'yticklabel', Months)

%     if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
    if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end
    
    if sp == 4
        ax = axes('visible', 'off');
        colormap(ax, cmaps_top)
        cbar = colorbar('horizontal');
        set(cbar,'xtick',cmaps_linspace_top, ...
            'XTickLabel',cmaps_y_label_top, ...
            'ticklength',cbar_tick_length, ...
            'XAxisLocation','bottom', 'fontsize',font_size)
        set(cbar,'units','centimeters','position', [...
            marg_l+x_sp*(cm(sp)-4)+gap_w*(cm(sp)-4), ...
            marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)-plot_cbar_gap, ...
            cbar_x*4+gap_w*3, ...
            cbar_y]);
        pointycbar(cbar)

        
    elseif sp == rowN*colN
        ax = axes('visible', 'off');
        colormap(ax, cmaps_bot);
        cbar = colorbar('horizontal');
        set(cbar,'xtick',cmaps_linspace_bot, ...
            'XTickLabel',cmaps_y_label_bot,...
            'ticklength',cbar_tick_length, ...
            'fontsize',font_size)
        set(cbar,'units','centimeters','position', [...
            (marg_l+x_sp*(cm(sp)-4)+gap_w*(cm(sp)-4)), ...
            (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)-plot_cbar_gap), ...
            cbar_x*4+gap_w*3, ...
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


%%
play(bird_f_path,[1 (get(bird_f_path, 'SampleRate')*3)]);

