%% fig 1: map of SBC currents
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_U_ek'])
load([data_path 'SACS_data/aus8_V_ek'])
load([data_path 'SACS_data/KDau_U_ek'])
load([data_path 'SACS_data/KDau_V_ek'])
load([data_path 'SACS_data/aus8_currents'])
load([data_path 'SACS_data/aus8_figures'])


%%
lon_u = aus8_coor.lon_u;
lon_v = aus8_coor.lon_v;
lat_u = aus8_coor.lat_u;
lat_v = aus8_coor.lat_v;

aus8_Uek = aus8_U_ek.mean;
KDau_Uek = KDau_U_ek.mean;
aus8_Vek = aus8_V_ek.mean;
KDau_Vek = KDau_V_ek.mean;

aus8_Vek_i2 = interp2(...
    lon_v, lat_v, aus8_Vek, lon_u, lat_u);
KDau_Vek_i2 = interp2(...
    lon_v, lat_v, KDau_Vek, lon_u, lat_u);

aus8_speed = sqrt(aus8_Uek.^2 + aus8_Vek_i2.^2);
KDau_speed = sqrt(KDau_Uek.^2 + KDau_Vek_i2.^2);


%%
screen_ratio = 0.75;
fig_n = 1;
rowcols = [2 1];
rowcols_size = [12 4.2]/screen_ratio; % cm
margs = [1.0 0.3 1.9 0.7]/screen_ratio; % cm
gaps = [0.4 0.7]/screen_ratio; % cm

magnif = 100;
cmap1_cont = [0 0.25 0.5 0.75 1 1.25 1.5 1.75 2]*magnif;
lvl_cmap1 = length(cmap1_cont)-1;
cmap1 = othercolor('Purples8', lvl_cmap1);
% cmap1(end,:) = [1 1 1];
cmaps_cont_length = length(cmap1_cont);
cmaps_linspace = linspace(0,1,cmaps_cont_length);
cmaps_y_label = cmap1_cont/magnif;

% lon_min = 115; lon_max = 147; lat_min = -47; lat_max = -32;
lon_min = 110; lon_max = 152; lat_min = -48; lat_max = -31;

x_chc = {aus8_coor.lon_u, aus8_coor.lon_u};
x_ind = [1 1];
y_chc = {aus8_coor.lat_u, aus8_coor.lat_u};

data = {aus8_speed*magnif, KDau_speed*magnif};

aus8_Uek_norm = aus8_Uek./(aus8_Uek.^2 + aus8_Vek_i2.^2).^(1/4);
aus8_Vek_norm = aus8_Vek_i2./(aus8_Uek.^2 + aus8_Vek_i2.^2).^(1/4);

KDau_Uek_norm = KDau_Uek./(KDau_Uek.^2 + KDau_Vek_i2.^2).^(1/4);
KDau_Vek_norm = KDau_Vek_i2./(KDau_Uek.^2 + KDau_Vek_i2.^2).^(1/4);

u_data = {aus8_Uek_norm*magnif, KDau_Uek_norm*magnif};
v_data = {aus8_Vek_norm*magnif, KDau_Vek_norm*magnif};

title_chc = ...
    {'ERA-Interim', ...
    'CORE2-NYF'};

for sp = 1 : rowcols(1)*rowcols(2)
    minmax{sp} = [cmap1_cont(1) cmap1_cont(end)];
    cmaps_custom{sp} = cmapcust(cmap1,cmap1_cont);
    
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
    
%     % create reference arrow
%     lat_ref = -34; lon_ref = 144;
%     lat_ref_ind = find(lat_u==lat_ref+1/16); lon_ref_ind = find(lon_u==lon_ref);
%     
%     data{sp}(lat_ref_ind-15:lat_ref_ind,lon_ref_ind:lon_ref_ind+32) = 0;
%     u_data{sp}(lat_ref_ind-15:lat_ref_ind,lon_ref_ind:lon_ref_ind+32) = 0;
%     v_data{sp}(lat_ref_ind-15:lat_ref_ind,lon_ref_ind:lon_ref_ind+32) = 0;
    
    if sp == 1
        s = 1.25;
        ref_magn = 1;
    elseif sp == 2
        s = 1;
        ref_magn = 1;
    end
    
%     data{sp}(lat_ref_ind-10 ,lon_ref_ind+6) = ref_magn*magnif;
%     u_data{sp}(lat_ref_ind-10 ,lon_ref_ind+6) = ref_magn*magnif;
    
    colormap(ax, cmaps_custom{sp});
    pcolor(x{sp}, y{sp}, data{sp})
    axis(axis_setup{sp})
    shading interp
    caxis([minmax{sp}(1) minmax{sp}(2)]);
    hold on
    
    for n = 1 : length(aus8_figures.cross.lon(1,:))
        plot(aus8_figures.cross.lon(1:2,n), ...
            aus8_figures.cross.lat(1:2,n), '--k', 'linewidth',1.2)
    end
    
    n = 7; m = 7;
    quiver(lon_u(1:n:end),lat_u(1:m:end),...
        u_data{sp}(1:m:end,1:n:end),...
        v_data{sp}(1:m:end,1:n:end),...
        s,'k')
    
%     text(lon_u(lon_ref_ind+6), lat_u(lat_ref_ind-9+5), ...
%         [num2str(ref_magn) ' $m^{2}/s$'], ...
%         'fontsize',font_size)
    
    h_tit = title(['(' lett(sp) ') Annual mean Ekman drift from ' ...
        title_chc{sp}], ...
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
        colormap(ax, cmap1);
        cbar = colorbar('horizontal');
        set(cbar,'xtick',cmaps_linspace, ...
            'XTickLabel',cmaps_y_label,...
            'fontsize',font_size,'ticklength',cbar_tick_length)
        set(cbar,'units','centimeters','position', [...
            marg_l, ...
            marg_b+plot_cbar_gap, ...
            cbar_x, ...
            cbar_y]);
        
        set(get(cbar,'xlabel'), ...
            'String','Ekman drift ($m^{2}/s$)', ...
            'fontsize',font_size)
        cbar.Label.Interpreter = 'latex';

%         pointycbar(cbar)
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
