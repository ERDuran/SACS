%% fig 2: meridional cross sections of zonal velocity differences
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_currents'])
load([data_path 'SACS_data/KDau_fcrt'])
Seasons = {'Summer', 'Autumn', 'Winter', 'Spring'};
load([data_path 'SACS_data/KDau_U_ek'])
load([data_path 'SACS_data/KDau_V_ek'])

%%
lon_u = aus8_coor.lon_u;
lon_v = aus8_coor.lon_v;
lat_u = aus8_coor.lat_u;
lat_v = aus8_coor.lat_v;

KDau_Uek_JFM = (KDau_U_ek.Jan + KDau_U_ek.Feb + KDau_U_ek.Mar)/3;
KDau_Vek_JFM = (KDau_V_ek.Jan + KDau_V_ek.Feb + KDau_V_ek.Mar)/3;
KDau_Vek_i2_JFM = interp2(...
    lon_v, lat_v, KDau_Vek_JFM, lon_u, lat_u);
KDau_speed_JFM = sqrt(KDau_Uek_JFM.^2 + KDau_Vek_i2_JFM.^2);

KDau_Uek_AMJ = (KDau_U_ek.Apr + KDau_U_ek.May + KDau_U_ek.Jun)/3;
KDau_Vek_AMJ = (KDau_V_ek.Apr + KDau_V_ek.May + KDau_V_ek.Jun)/3;
KDau_Vek_i2_AMJ = interp2(...
    lon_v, lat_v, KDau_Vek_AMJ, lon_u, lat_u);
KDau_speed_AMJ = sqrt(KDau_Uek_AMJ.^2 + KDau_Vek_i2_AMJ.^2);

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

LCE_u = lon_u >= 115 & lon_u <= 124;
zSAC_u = lon_u >= 124 & lon_u <= 132;
sSAC_u = lon_u >= 132 & lon_u <= 141;
ZC_u = lon_u >= 141 & lon_u <= 147;

SBC_Ut_mean(:,1) = mean(SBC_Ut(:,LCE_u),2);
SBC_Ut_mean(:,2) = mean(SBC_Ut(:,zSAC_u),2);
SBC_Ut_mean(:,3) = mean(SBC_Ut(:,sSAC_u),2);
SBC_Ut_mean(:,4) = mean(SBC_Ut(:,ZC_u),2)

DRC_Ut_mean(:,1) = mean(DRC_Ut(:,LCE_u),2);
DRC_Ut_mean(:,2) = mean(DRC_Ut(:,zSAC_u),2);
DRC_Ut_mean(:,3) = mean(DRC_Ut(:,sSAC_u),2);
DRC_Ut_mean(:,4) = mean(DRC_Ut(:,lon_u >= 141 & lon_u <= 144),2)


%%
screen_ratio = 0.75;
fig_n = 1;
rowcols = [2 2];
rowcols_size = [14 6]/screen_ratio/2; % cm
margs = [1 0.2 1.8 0.6]/screen_ratio; % cm
gaps = [0.7 0.7]/screen_ratio; % cm
plot_cbar_gap = 0.7/screen_ratio;
cbar_x = rowcols_size(1);
cbar_y = 0.2/screen_ratio;

magnif = 100;
cmap1_cont = [0 0.25 0.5 0.75 1 1.25 1.5 1.75 2]*magnif;
lvl_cmap1 = length(cmap1_cont)-1;
cmap1 = othercolor('Purples8', lvl_cmap1);
% cmap1(end,:) = [1 1 1];
cmaps_cont_length = length(cmap1_cont);
cmaps_linspace = linspace(0,1,cmaps_cont_length);
cmaps_y_label = cmap1_cont/magnif;

lon_min = 115; lon_max = 147;
lat_min = {-1.5,-20,-48,-48}; lat_max = {2,10,-31,-31};

x_chc = {aus8_coor.lon_u, aus8_coor.lon_v};
x_ind = [1 1 1 1];
y_chc = {aus8_coor.lat_u, aus8_coor.lat_u};

data = {SBC_Ut, DRC_Ut, KDau_speed_JFM*magnif, KDau_speed_AMJ*magnif};
data_mean = {KDau_fcrt.SBC_Ut.mean, KDau_fcrt.DRC_Ut.mean, ...
    KDau_fcrt.SBC_Wtc.mean};

KDau_Uek_JFM_norm = KDau_Uek_JFM./(KDau_Uek_JFM.^2 + KDau_Vek_i2_JFM.^2).^(1/4);
KDau_Vek_i2_JFM_norm = KDau_Vek_i2_JFM./(KDau_Uek_JFM.^2 + KDau_Vek_i2_JFM.^2).^(1/4);

KDau_Uek_AMJ_norm = KDau_Uek_AMJ./(KDau_Uek_AMJ.^2 + KDau_Vek_i2_AMJ.^2).^(1/4);
KDau_Vek_i2_AMJ_norm = KDau_Vek_i2_AMJ./(KDau_Uek_AMJ.^2 + KDau_Vek_i2_AMJ.^2).^(1/4);

u_data = {KDau_Uek_JFM_norm*magnif, KDau_Uek_AMJ_norm*magnif};
v_data = {KDau_Vek_i2_JFM_norm*magnif, KDau_Vek_i2_AMJ_norm*magnif};

data_label = ...
    {'$\mathcal{U}_{SBC}$', '$\mathcal{U}_{FC}$', '$\mathcal{W}_{SBC}$'};

% title_chc = {'U''', 'V''', 'U_{g}''', 'V_{g}'''};

for sp = 1 : rowcols(1)*rowcols(2)
    axis_setup{sp} = [lon_min lon_max lat_min{sp} lat_max{sp}];
    x{sp} = x_chc{x_ind(sp)};
    
    minmax{sp} = [cmap1_cont(1) cmap1_cont(end)];
    cmaps_custom{sp} = cmapcust(cmap1,cmap1_cont);
    
    y{sp} = y_chc{x_ind(sp)};
    
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

desired_length = 0.05/screen_ratio; %cm
if y_sp > x_sp, long_side = y_sp; else, long_side = x_sp; end
norm_length = desired_length/long_side;
fig_tick_length = [norm_length; 0.01];
if cbar_y > cbar_x, long_side = cbar_y; else, long_side = cbar_x; end
norm_length = desired_length/long_side;
cbar_tick_length = norm_length;

line_width = 0.5;
title_chc = ...
    {'Summer', ...
    'Autumn'};

for sp = 1 : rowN*colN
    subplot_x = marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1);
    subplot_y = marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1);
    ax = axes('Units','centimeters', ...
        'Position',[subplot_x,subplot_y,x_sp,y_sp]);
    
     if sp <= 2   
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
        
        if sp == 1
            text(119, 1.8, num2str(round(SBC_Ut_mean(2,1)*10)/10), 'color', hhc4, ...
                'fontsize',font_size)
            text(127, 1.8, num2str(round(SBC_Ut_mean(2,2)*10)/10), 'color', hhc4, ...
                'fontsize',font_size)
            text(135, 1.8, num2str(round(SBC_Ut_mean(2,3)*10)/10), 'color', hhc4, ...
                'fontsize',font_size)
            text(143, 1.8, num2str(round(SBC_Ut_mean(2,4)*10)/10), 'color', hhc4, ...
                'fontsize',font_size)
            
            text(119, -1.3, num2str(round(SBC_Ut_mean(1,1)*10)/10), 'color', hhc2, ...
                'fontsize',font_size)
            text(127, -1.3, num2str(round(SBC_Ut_mean(1,2)*10)/10), 'color', hhc2, ...
                'fontsize',font_size)
            text(135, -1.3, num2str(round(SBC_Ut_mean(1,3)*10)/10), 'color', hhc2, ...
                'fontsize',font_size)
            text(143, -1.3, num2str(round(SBC_Ut_mean(1,4)*10)/10), 'color', hhc2, ...
                'fontsize',font_size)
            
        elseif sp == 2
            text(119, 7.5, num2str(round(DRC_Ut_mean(2,1)*10)/10), 'color', hhc4, ...
                'fontsize',font_size)
            text(127, 7.5, num2str(round(DRC_Ut_mean(2,2)*10)/10), 'color', hhc4, ...
                'fontsize',font_size)
            text(135, 7.5, num2str(round(DRC_Ut_mean(2,3)*10)/10), 'color', hhc4, ...
                'fontsize',font_size)
            text(143, 7.5, num2str(round(DRC_Ut_mean(2,4)*10)/10), 'color', hhc4, ...
                'fontsize',font_size)
            
            text(119, -17.5, num2str(round(DRC_Ut_mean(1,1)*10)/10), 'color', hhc2, ...
                'fontsize',font_size)
            text(127, -17.5, num2str(round(DRC_Ut_mean(1,2)*10)/10), 'color', hhc2, ...
                'fontsize',font_size)
            text(135, -8.5, num2str(round(DRC_Ut_mean(1,3)*10)/10), 'color', hhc2, ...
                'fontsize',font_size)
            text(143, -8.5, num2str(round(DRC_Ut_mean(1,4)*10)/10), 'color', hhc2, ...
                'fontsize',font_size)
        end
        
        if sp == 2
            %         MTH_l = Seasons;
            MTH_l = {'Summer', 'Autumn', 'Winter', 'Spring', 'Annual mean'};
            h_leg = legend(MTH_l);
            set(h_leg,'units','centimeters', 'orientation','vertical', ...
                'fontsize',font_size)
            h_pos1 = get(h_leg,'position')/screen_ratio;
            set(h_leg,'position', [...
                (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1))+5, ...
                (marg_b+y_sp*(cm(sp)-1)+gap_h*(cm(sp)-1))-0.75, ...
                h_pos1(3), ...
                h_pos1(4)]);
        end
        
        if sp == 2
            set(ax,'layer','top','color',nan_color,...
                'ticklength',fig_tick_length, ...
                'fontsize',font_size,'tickdir','out', ...
                'xtick', x_tick, 'xticklabel', x_tick_label, ...
                'ytick', -20:5:10)
        elseif sp == 3
            set(ax,'layer','top','color',nan_color,...
                'ticklength',fig_tick_length, ...
                'fontsize',font_size,'tickdir','out', ...
                'xtick', x_tick, 'xticklabel', x_tick_label, ...
                'ytick', lat_min{sp}:0.5:lat_max{sp})
        else
            set(ax,'layer','top','color',nan_color,...
                'ticklength',fig_tick_length, ...
                'fontsize',font_size,'tickdir','out', ...
                'xtick', x_tick, 'xticklabel', x_tick_label, ...
                'ytick', -2:0.5:2)
        end
        
        if row_ind(sp) ~= rowN
            set(gca,'xticklabel','')
        else
            xlabel('Longitude ($^{\circ}E$)')
        end
%         if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end
        if sp == 1
            ylabel('Transport ($Sv$)')
        end
        
    else
        % create reference arrow
%         lat_ref = -34; lon_ref = 142;
%         lat_ref_ind = find(lat_u==lat_ref+1/16); lon_ref_ind = find(lon_u==lon_ref);
%         
%         data{sp}(lat_ref_ind-15:lat_ref_ind+2,lon_ref_ind:lon_ref_ind+32) = 0;
%         u_data{sp-2}(lat_ref_ind-15:lat_ref_ind,lon_ref_ind:lon_ref_ind+32) = 0;
%         v_data{sp-2}(lat_ref_ind-15:lat_ref_ind,lon_ref_ind:lon_ref_ind+32) = 0;
        
        if sp == 3
            s = 1.15;
            ref_magn = 1;
        elseif sp == 4
            s = 0.75;
            ref_magn = 1;
        end
        
%         data{sp}(lat_ref_ind-10 ,lon_ref_ind+8) = ref_magn*magnif;
%         u_data{sp-2}(lat_ref_ind-11 ,lon_ref_ind+8) = ref_magn*magnif;
        
        colormap(ax, cmaps_custom{sp});
        pcolor(x{sp}, y{sp}, data{sp})
        axis(axis_setup{sp})
        shading interp
        caxis([minmax{sp}(1) minmax{sp}(2)]);
        hold on
        
        %     for n = 1 : length(aus8_figures.cross.lon(1,:))
        %         plot(aus8_figures.cross.lon(1:2,n), ...
        %             aus8_figures.cross.lat(1:2,n), '--k', 'linewidth',1.2)
        %     end
        
        n = 10; m = 10;
        quiver(lon_u(1:n:end),lat_u(1:m:end),...
            u_data{sp-2}(1:m:end,1:n:end),...
            v_data{sp-2}(1:m:end,1:n:end),...
            s,'k')
        
%         text(lon_u(lon_ref_ind+6), lat_u(lat_ref_ind-9+5), ...
%             [num2str(ref_magn) ' $m^{2}/s$'], ...
%             'fontsize',font_size)
        
        h_tit = title(['(' lett(sp) ') ' title_chc{sp-2} ...
            ' Ekman drift in CORE2-NYF'], ...
            'horizontalalignment','left', 'fontsize',font_size, ...
            'Interpreter','latex');
        h_tit.Position(1) = axis_setup{sp}(1);
        grid
        nan_color = [0.7 0.7 0.7];
        set(ax,'layer','top','color',nan_color,...
            'fontsize',font_size,'tickdir','out', ...
            'ticklength',fig_tick_length, ...
            'xtick', lon_min:2:lon_max, 'ytick', lat_min{3}:2:lat_max{3})
        if row_ind(sp) ~= rowN
            set(gca,'xticklabel','')
            
        else
            xlabel('Longitude ($^{\circ}$E)', 'fontsize', font_size)
        end
        
        if sp == 3
            ylabel('Latitude ($^{\circ}$N)', 'fontsize', font_size)
        end
        if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end
        
        
        
        if sp == rowN*colN
            ax = axes('visible', 'off');
            colormap(ax, cmap1);
            cbar = colorbar('horizontal');
            set(cbar,'xtick',cmaps_linspace, ...
                'XTickLabel',cmaps_y_label,...
                'fontsize',font_size,'ticklength',cbar_tick_length)
            set(cbar,'units','centimeters','position', [...
                marg_l, ...
                marg_b+plot_cbar_gap-2.25, ...
                cbar_x*2+gap_w, ...
                cbar_y]);
            
            set(get(cbar,'xlabel'), ...
                'String','Ekman drift ($m^{2}/s$)', ...
                'fontsize',font_size)
            cbar.Label.Interpreter = 'latex';
            
            %         pointycbar(cbar)
        end
        
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
    '-dpng', '-r400')
% print(fig, ...
%     ['~/Duran2017/SACS/10319442jbhpxfsdfvwy/' scriptname(1:3) ...
%     '_fig' num2str(fig_n) '_'], ...
%     '-dpng', '-r400')


%%
play(bird_f_path,[1 (get(bird_f_path, 'SampleRate')*3)]);

