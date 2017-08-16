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
% Months{13} = 'mean';

SBC_Ut = NaN(length(Months), length(lon_u));
[SBC_Vtnc, SBC_Vtsc, SBC_Wtc_real] = ...
    deal(NaN(length(Months), length(lon_v)));

for t = 1 : 12
    SBC_Ut(t,:) = aus8_currents.SBC_Ut.(Months{t});
    SBC_Vtnc(t,:) = aus8_currents.SBC_Vtnc.(Months{t});
    SBC_Vtsc(t,:) = aus8_currents.SBC_Vtsc.(Months{t});
    SBC_Wtc_real(t,:) = aus8_currents.SBC_Wtc_real.(Months{t});
    
end

%%
screen_ratio = 0.75;
fig_n = 1;
rowcols = [1 4];
rowcols_size = [4 4]/screen_ratio; % cm
margs = [0.6 1.2 0.6 0.6]/screen_ratio; % cm
gaps = [0.2 0.8]/screen_ratio; % cm
plot_cbar_gap = 0.3/screen_ratio;
cbar_x = rowcols_size(1);
cbar_y = 0.2/screen_ratio;

magnif = 100;
cmap1_cont = -[3 2 1.5 1 0.75 0.5 0.25 0.1 0]*magnif;
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

lon_min = 115; lon_max = 147; lat_min = 1; lat_max = 12;

x_chc = {aus8_coor.lon_u, aus8_coor.lon_v};
x_ind = [1 2 2 2];

data = {SBC_Ut*magnif, SBC_Vtnc*magnif, ...
    SBC_Vtsc*magnif, SBC_Wtc_real*magnif};

title_chc = {'U''', 'V''', 'U_{g}''', 'V_{g}'''};

for sp = 1 : rowcols(1)*rowcols(2)
    minmax{sp} = [cmaps_cont(1) cmaps_cont(end)];
    cmaps_custom{sp} = cmapcust(cmaps,cmaps_cont);
    
    axis_setup{sp} = [lon_min lon_max lat_min lat_max];
    x{sp} = x_chc{x_ind(sp)};
    y{sp} = 1:12;

end

lett = 'a':'z';
font_size = 8*screen_ratio;
fig_color = [0.7 0.7 0.7];

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
set(fig,'units','centimeters','paperunits','centimeters', ...
    'inverthardcopy','off','color',[1 1 1],...
    'paperposition',[0 0 ...
    (marg_l+colN*x_sp+gap_w*(colN-1)+marg_r) ...
    (marg_b+rowN*y_sp+gap_h*(rowN-1)+marg_t)]*screen_ratio,...
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
    
    h_tit = title(['(']);
        h_tit.Position(1) = axis_setup{sp}(1);
    grid
    x_tick = lon_min:2:lon_max;
    x_tick_label = num2cell(x_tick);
    for c = 1:2:17
        x_tick_label{c} = '';
    end
    
    set(ax,'layer','top','color',fig_color,...
        'fontsize',font_size,'tickdir','out', ...
        'xtick', x_tick, 'xticklabel', x_tick_label, ...
        'ytick', lat_min:1:lat_max, 'yticklabel', Months)
    

    if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
    if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end
    
    if sp == rowN*colN
        ax = axes;
        set(gca,'Visible','off')
        colormap(ax, cmaps);
        cbar = colorbar('horizontal');
        set(cbar,'xtick',cmaps_linspace, ...
            'XTickLabel',cmaps_y_label,...
            'fontsize',font_size)
        set(cbar,'units','centimeters','position', [...
            (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)+plot_cbar_gap), ...
            (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
            cbar_x, ...
            cbar_y]);
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
close all


%%
play(bird_f_path,[1 (get(bird_f_path, 'SampleRate')*3)]);

