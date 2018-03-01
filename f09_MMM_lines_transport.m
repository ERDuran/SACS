%% fig 2: meridional cross sections of zonal velocity differences
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_currents'])
load([data_path 'SACS_data/KDau_fcrt'])
Seasons = {'Summer', 'Autumn', 'Winter', 'Spring'};


%%
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;
MTH = aus8_coor.MTH;

SBC_Ut = NaN(length(MTH), length(lon_u));
SBC_Wtc = ...
    NaN(length(MTH), length(lon_v));
DRC_Ut = NaN(length(MTH), length(lon_u));
for t = 1 : 4
%     SBC_Ut(t,:) = aus8_currents.SBC_Ut.(MTH{t});
%     SBC_Vtnc(t,:) = aus8_currents.SBC_Vtnc.(MTH{t});
%     SBC_Vtsc(t,:) = aus8_currents.SBC_Vtsc.(MTH{t});
%     SBC_Wtc_real(t,:) = aus8_currents.SBC_Wtc_real.(MTH{t});
%     DRC_Ut(t,:) = aus8_currents.DRC_Ut.(MTH{t});
%     DRC_Vtnc(t,:) = aus8_currents.DRC_Vtnc.(MTH{t});
%     DRC_Vtsc(t,:) = aus8_currents.DRC_Vtsc.(MTH{t});
%     DRC_Wtc_real(t,:) = aus8_currents.DRC_Wtc_real.(MTH{t});
    
    SBC_Ut(t,:) = KDau_fcrt.SBC_Ut.(MTH{t});
    SBC_Wtc(t,:) = KDau_fcrt.SBC_Wtc.(MTH{t});
    DRC_Ut(t,:) = KDau_fcrt.DRC_Ut.(MTH{t});
end


%%
screen_ratio = 0.75;
fig_n = 1;
rowcols = [3 1];
rowcols_size = [14 6]/screen_ratio/2; % cm
margs = [1 0.2 0.8 0.6]/screen_ratio; % cm
gaps = [0.4 0.8]/screen_ratio; % cm
plot_cbar_gap = 0.7/screen_ratio;
cbar_x = rowcols_size(1);
cbar_y = 0.2/screen_ratio;


lon_min = 115; lon_max = 147; 
lat_min = {-1.5,-20,-4}; lat_max = {2,7.5,0};

x_chc = {aus8_coor.lon_u, aus8_coor.lon_v};
x_ind = [1 1 2];

data = {SBC_Ut, DRC_Ut, SBC_Wtc};
data_mean = {KDau_fcrt.SBC_Ut.mean, KDau_fcrt.DRC_Ut.mean, ...
    KDau_fcrt.SBC_Wtc.mean};

data_label = ...
    {'$\mathcal{U}_{SBCs}$', '$\mathcal{U}_{FC}$', '$\mathcal{W}_{SBCs}$'};

% title_chc = {'U''', 'V''', 'U_{g}''', 'V_{g}'''};

for sp = 1 : rowcols(1)*rowcols(2)
    axis_setup{sp} = [lon_min lon_max lat_min{sp} lat_max{sp}];
    x{sp} = x_chc{x_ind(sp)};

end

lett = 'a':'z';
font_size = 8*screen_ratio;
nan_color = [1 1 1];
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

set(fig,'units','centimeters','paperunits','centimeters', ...
    'inverthardcopy','off','color',fig_color,...
    'paperposition',[0 0 fig_x fig_y]*screen_ratio,...
    'position',[0 0 fig_x fig_y]);

line_width = 0.5;

for sp = 1 : rowN*colN
    subplot_x = marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1);
    subplot_y = marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1);
    ax = axes('Units','centimeters', ...
        'Position',[subplot_x,subplot_y,x_sp,y_sp]);
    
    hh = plot([124 124], [-30 30],...
        'linewidth', line_width+0.2);
    hhc1 = get(hh,'color');
    set(hh, 'color', [0.7 0.7 0.7])
    set(get(get(hh,'Annotation'),'LegendInformation'), ...
        'IconDisplayStyle','off');
    hold on
    
    hh = plot([132 132], [-30 30], 'linestyle', '--', ...
        'linewidth', line_width+0.2);
    hhc2 = get(hh,'color');
    set(hh, 'color', [0.7 0.7 0.7])
    set(get(get(hh,'Annotation'),'LegendInformation'), ...
        'IconDisplayStyle','off');
    
    hh = plot([141 141], [-30 30], ...
        'linewidth', line_width+0.2);
    hhc3 = get(hh,'color');
    set(hh, 'color', [0.7 0.7 0.7])
    set(get(get(hh,'Annotation'),'LegendInformation'), ...
        'IconDisplayStyle','off');
    
    hh = plot([1 1], [1 1]);
    hhc4 = get(hh,'color');
    set(get(get(hh,'Annotation'),'LegendInformation'), ...
        'IconDisplayStyle','off');
    
    plot(x{sp}(1,:), data{sp}(1,:), 'color', hhc2)
    plot(x{sp}(1,:), data{sp}(2,:), 'color', hhc4)
    plot(x{sp}(1,:), data{sp}(3,:), 'color', hhc1)
    plot(x{sp}(1,:), data{sp}(4,:), 'color', hhc3)
    
    plot(x{sp}, data_mean{sp}, 'k--')
    
    plot([124 124], [-30 30], ...
    'color', [0.5 0.5 0.5], 'linewidth', line_width)
    plot([141 141], [-30 30], ...
    'color', [0.5 0.5 0.5], 'linewidth', line_width)

    axis(axis_setup{sp})
    
    h_tit = title(['(' lett(sp) ') ' data_label{sp}], ...
        'horizontalalignment','left', 'fontsize',font_size);
    h_tit.Position(1) = axis_setup{sp}(1);
    grid
    x_tick = lon_min:2:lon_max;
    x_tick_label = num2cell(x_tick);

    if sp == 3
        MTH_l = Seasons;
        MTH_l{5} = 'mean';
        h_leg = legend(MTH_l);
        set(h_leg,'units','centimeters', 'orientation','vertical', ...
            'fontsize',font_size)
        h_pos1 = get(h_leg,'position')/screen_ratio;
        set(h_leg,'position', [...
            (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)), ...
            (marg_b+y_sp*(cm(sp)-1)+gap_h*(cm(sp)-1)), ...
            h_pos1(3), ...
            h_pos1(4)]);
    end
    
    if sp == 2
        set(ax,'layer','top','color',nan_color,...
            'ticklength',fig_tick_length, ...
            'fontsize',font_size,'tickdir','out', ...
            'xtick', x_tick, 'xticklabel', x_tick_label, ...
            'ytick', lat_min{sp}:2.5:lat_max{sp})
    else
        set(ax,'layer','top','color',nan_color,...
            'ticklength',fig_tick_length, ...
            'fontsize',font_size,'tickdir','out', ...
            'xtick', x_tick, 'xticklabel', x_tick_label, ...
            'ytick', lat_min{sp}:0.5:lat_max{sp})
    end

    if row_ind(sp) ~= rowN
        set(gca,'xticklabel','')
    else
        xlabel('Longitude')
    end
    if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end
    ylabel('Transport ($Sv$)')
    
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


%%
play(bird_f_path,[1 (get(bird_f_path, 'SampleRate')*3)]);

