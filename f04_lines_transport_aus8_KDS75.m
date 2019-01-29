%% fig 2: meridional cross sections of zonal velocity differences
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_currents'])
load([data_path 'SACS_data/KDau_currents'])
load([data_path 'SACS_data/KDau_fcrt'])

lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;


%%
close all
screen_ratio = 0.75;
fig_n = 1;
rowcols = [2 2];
rowcols_size = [7 4.8]/screen_ratio; % cm
margs = [1 0.2 2.0 0.4]/screen_ratio; % cm
gaps = [0.8 0.6]/screen_ratio; % cm
plot_cbar_gap = 0.5/screen_ratio;
cbar_x = 0.3/screen_ratio;
cbar_y = rowcols_size(2);
lett = 'a':'z';

font_size = 8*screen_ratio;
fig_color = [1 1 1];

fig = figure(fig_n);
rowN = rowcols(1); colN = rowcols(2);
[rm, cm] = meshgrid(rowN:-1:1, 1:colN);
x_sp = rowcols_size(1); % x subplot length
y_sp = rowcols_size(2); % y subplot length
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = gaps(1); % gap width between subplots
gap_h = gaps(2); % gap height between subplots
marg_b = margs(3); % bottom_aus8 margin
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

line_width = 0.5;
dashed = 0.1;
big_line = 0;

sp = 1;
subplot_x = marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1);
subplot_y = marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1);
ax = axes('Units','centimeters', ...
    'Position',[subplot_x,subplot_y,x_sp,y_sp]);

hh = plot([124 124], [-30 30],...
    'color', [0.7 0.7 0.7], 'linewidth', line_width+0.2);
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');
hold on
hh = plot([132 132], [-30 30], 'linestyle', '--', ...
    'color', [0.7 0.7 0.7], 'linewidth', line_width+0.2);
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');
hh = plot([141 141], [-30 30], ...
    'color', [0.8 0.8 0.8], 'linewidth', line_width+0.2);
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');

plot(lon_u, aus8_currents.SBC_Ut.mean, 'k-', ...
    'linewidth', line_width+big_line)
plot(lon_u, KDau_currents.SBC_Ut.mean, 'k--', ...
    'linewidth', line_width-dashed)
plot(lon_u, KDau_fcrt.SBC_Ut.mean, 'k:', 'linewidth', line_width)

plot(lon_v, aus8_currents.SBC_Vtnc.mean, 'g-', ...
    'linewidth', line_width+big_line)
plot(lon_v, KDau_currents.SBC_Vtnc.mean, 'g--', ...
    'linewidth', line_width-dashed)
plot(lon_v, KDau_fcrt.SBC_Vtnc.mean, 'g:', 'linewidth', line_width)

text(119, 1.9, 'LCE', 'fontsize',font_size)
text(131.1, 1.9, 'SAC', 'fontsize',font_size)
text(143, 1.9, 'ZC', 'fontsize',font_size)
text(131, 1.6, 'Zonal shelf', 'fontsize',font_size, ...
    'horizontalalignment', 'right')
text(133, 1.6, 'Slanted shelf', 'fontsize',font_size)

axis([115 147 -2.1 2.1])
h_tit = title(['(' lett(sp) ...
    ') SBC long-shore and CC cumulative transports'], ...
    'fontsize',font_size, 'horizontalalignment','left');
h_tit.Position(1) = 115;

grid
h_leg = legend(...
    'ZD CARS $\mathcal{U}_{SBC}$', ...
    'ZD MOM01 $\mathcal{U}_{SBC}$', ...
    'Full MOM01 $\mathcal{U}_{SBC}$', ...
    'ZD CARS $\mathcal{V}_{CC}$', ...
    'ZD MOM01 $\mathcal{V}_{CC}$', ...
    'Full MOM01 $\mathcal{V}_{CC}$');

set(h_leg,'units','centimeters', 'orientation','vertical', ...
    'fontsize',font_size)
title(h_leg,['Legend (' lett(sp) ')'],'fontsize',font_size)
h_pos1 = get(h_leg,'position')/screen_ratio;
set(h_leg,'position', [...
    (marg_l+x_sp*(cm(sp))+gap_w*(cm(sp))), ...
    (marg_b+y_sp*(cm(sp)+1)+gap_h*(cm(sp))-h_pos1(4)), ...
    h_pos1(3), ...
    h_pos1(4)]);

set(ax,'layer','top','color',fig_color,...
    'fontsize',font_size,'tickdir','out',...
    'ticklength',fig_tick_length, ...
    'xtick',115:2:155,'ytick',-3:0.5:3)
if row_ind(sp) ~= rowN
    set(gca,'xticklabel','')
else
    xlabel('Longitude')
end
if col_ind(sp) == 1
    ylabel('Transport ($Sv$)')
end


sp = 3;
subplot_x = marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1);
subplot_y = marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1);
ax = axes('Units','centimeters', ...
    'Position',[subplot_x,subplot_y,x_sp,y_sp]);

hh = plot([124 124], [-30 30],...
    'color', [0.7 0.7 0.7], 'linewidth', line_width+0.2);
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');
hold on
hh = plot([132 132], [-30 30], 'linestyle', '--', ...
    'color', [0.7 0.7 0.7], 'linewidth', line_width+0.2);
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');
hh = plot([141 141], [-30 30], ...
    'color', [0.8 0.8 0.8], 'linewidth', line_width+0.2);
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');

plot(lon_v, aus8_currents.SBC_Vtsc.mean, 'r-', ...
    'linewidth', line_width+big_line)
plot(lon_v, KDau_currents.SBC_Vtsc.mean, 'r--', ...
    'linewidth', line_width-dashed)
plot(lon_v, KDau_fcrt.SBC_Vtsc.mean, 'r:', 'linewidth', line_width)

plot(lon_v, aus8_currents.SBC_Wtc_real.mean, 'b-', ...
    'linewidth', line_width+big_line)
plot(lon_v, KDau_currents.SBC_Wtc_real.mean, 'b--', ...
    'linewidth', line_width-dashed)
plot(lon_v, KDau_fcrt.SBC_Wtc.mean, 'b:', 'linewidth', line_width)

axis([115 147 -2.1 2.1])
h_tit = title(['(' lett(sp-1) ...
    ') SBC cumulative transports with the FC'], ...
    'fontsize',font_size, 'horizontalalignment','left');
h_tit.Position(1) = 115;

grid
h_leg = legend(...
    'ZD CARS $\mathcal{V}_{SBC}$', ...
    'ZD MOM01 $\mathcal{V}_{SBC}$', ...
    'Full MOM01 $\mathcal{V}_{SBC}$', ...
    'ZD CARS $\mathcal{W}_{SBC}$', ...
    'ZD MOM01 $\mathcal{W}_{SBC}$', ...
    'Full MOM01 $\mathcal{W}_{SBC}$');

set(h_leg,'units','centimeters', 'orientation','vertical', ...
    'fontsize',font_size)
title(h_leg,['Legend (' lett(sp-1) ')'],'fontsize',font_size)
h_pos1 = get(h_leg,'position')/screen_ratio;
set(h_leg,'position', [...
    (marg_l+x_sp*(cm(sp))+gap_w*(cm(sp))), ...
    (marg_b+y_sp*(cm(sp))+gap_h*(cm(sp))), ...
    h_pos1(3), ...
    h_pos1(4)]);

set(ax,'layer','top','color',fig_color,...
    'fontsize',font_size,'tickdir','out',...
    'ticklength',fig_tick_length, ...
    'xtick',115:2:155,'ytick',-3:0.5:3)
if col_ind(sp) == 1
    ylabel('Transport ($Sv$)')
end


small_arr = 3;
% Cape Leeuwin
hh = arrow([115,-2.9], [115,-2.45], small_arr, ...
    'Facecolor', 'k', 'edgecolor', 'k');
text(115,-2.9,'Cape Leeuwin', 'fontsize', font_size, ...
    'color', 'k','Rotation',-45)
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');
% Albany
hh = arrow([118,-2.9], [118,-2.45], small_arr, ...
    'Facecolor', 'k', 'edgecolor', 'k');
text(118,-2.9,'Albany', 'fontsize', font_size, ...
    'color', 'k','Rotation',-45)
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');
% Cape Pasley
hh = arrow([123.4,-2.9], [123.4,-2.45], small_arr, ...
    'Facecolor', 'k', 'edgecolor', 'k');
text(123.4,-2.9,'Cape Pasley', 'fontsize', font_size, ...
    'color', 'k','Rotation',-45)
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');
% West GAB
hh = arrow([128,-2.9], [128,-2.45], small_arr, ...
    'Facecolor', 'k', 'edgecolor', 'k');
text(128,-2.9,'West GAB', 'fontsize', font_size, ...
    'color', 'k','Rotation',-45)
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');  
% East GAB
hh = arrow([132,-2.9], [132,-2.45], small_arr, ...
    'Facecolor', 'k', 'edgecolor', 'k');
text(132,-2.9,'East GAB', 'fontsize', font_size, ...
    'color', 'k','Rotation',-45)
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');  
% Cape Carnot
hh = arrow([135.7,-2.9], [135.7,-2.45], small_arr, ...
    'Facecolor', 'k', 'edgecolor', 'k');
text(135.7,-2.9,'Cape Carnot', 'fontsize', font_size, ...
    'color', 'k','Rotation',-45)
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');  
% Portland
hh = arrow([141.3,-2.9], [141.3,-2.45], small_arr, ...
    'Facecolor', 'k', 'edgecolor', 'k');
text(141.3,-2.9,'Portland', 'fontsize', font_size, ...
    'color', 'k','Rotation',-45)
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');  
%  King Is.
hh = arrow([144,-2.9], [144,-2.45], small_arr, ...
    'Facecolor', 'k', 'edgecolor', 'k');
text(144,-2.9,' King Island', 'fontsize', font_size, ...
    'color', 'k','Rotation',-45)
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');  
% South East Cape
hh = arrow([146.7,-2.9], [146.7,-2.45], small_arr, ...
    'Facecolor', 'k', 'edgecolor', 'k');
text(146.7,-2.9,'South East C.', 'fontsize', font_size, ...
    'color', 'k','Rotation',-45)
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');  

sp = 4;
subplot_x = marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1);
subplot_y = marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1);
ax = axes('Units','centimeters', ...
    'Position',[subplot_x,subplot_y,x_sp,y_sp]);

hh = plot([124 124], [-30 30],...
    'color', [0.7 0.7 0.7], 'linewidth', line_width+0.2);
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');
hold on
hh = plot([132 132], [-30 30], 'linestyle', '--', ...
    'color', [0.7 0.7 0.7], 'linewidth', line_width+0.2);
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');
hh = plot([141 141], [-30 30], ...
    'color', [0.8 0.8 0.8], 'linewidth', line_width+0.2);
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');

plot(lon_u, aus8_currents.DRC_Ut.mean, 'k-', ...
    'linewidth', line_width+big_line)
plot(lon_u, KDau_currents.DRC_Ut.mean, 'k--', ...
    'linewidth', line_width-dashed)
plot(lon_u, KDau_fcrt.DRC_Ut.mean, 'k:', 'linewidth', line_width)
plot(lon_v, aus8_currents.DRC_Vtnc.mean, 'r-', ...
    'linewidth', line_width+big_line)
plot(lon_v, KDau_currents.DRC_Vtnc.mean, 'r--', ...
    'linewidth', line_width-dashed)
plot(lon_v, KDau_fcrt.DRC_Vtnc.mean, 'r:', 'linewidth', line_width)
plot(lon_v, aus8_currents.DRC_Vtsc.mean, 'm-', ...
    'linewidth', line_width+big_line)
plot(lon_v, KDau_currents.DRC_Vtsc.mean, 'm--', ...
    'linewidth', line_width-dashed)
plot(lon_v, KDau_fcrt.DRC_Vtsc.mean, 'm:', 'linewidth', line_width)
plot(lon_v, aus8_currents.DRC_Wtc_real.mean, 'b-', ...
    'linewidth', line_width+big_line)
plot(lon_v, KDau_currents.DRC_Wtc_real.mean, 'b--', ...
    'linewidth', line_width-dashed)
plot(lon_v, KDau_fcrt.DRC_Wtc.mean, 'b:', 'linewidth', line_width)
plot(lon_v, KDau_fcrt.BOT_Wtc.mean, 'c:', 'linewidth', line_width)

axis([115 147 -16.2 2.6])
h_tit = title(['(' lett(sp-1) ...
    ') FC long-shore and cumulative horiz. and vert. transports'], ...
    'fontsize',font_size, 'horizontalalignment','left');
h_tit.Position(1) = 115;

grid
h_leg = legend(...
    'ZD CARS $\mathcal{U}_{FC}$', ...
    'ZD MOM01 $\mathcal{U}_{FC}$', ...
    'Full MOM01 $\mathcal{U}_{FC}$', ...
    'ZD CARS $\mathcal{V}_{FC}$', ...
    'ZD MOM01 $\mathcal{V}_{FC}$', ...
    'Full MOM01 $\mathcal{V}_{FC}$', ...
    'ZD CARS $\mathcal{V}_{OF}$', ...
    'ZD MOM01 $\mathcal{V}_{OF}$', ...
    'Full MOM01 $\mathcal{V}_{OF}$', ...
    'ZD CARS $\mathcal{W}_{FC}$', ...
    'ZD MOM01 $\mathcal{W}_{FC}$', ...
    'Full MOM01 $\mathcal{W}_{FC}$', ...
    'Full MOM01 $\mathcal{W}_{bot}$');
set(h_leg,'units','centimeters', 'orientation','vertical', ...
    'fontsize',font_size)
title(h_leg,['Legend (' lett(sp-1) ')'],'fontsize',font_size)
h_pos1 = get(h_leg,'position')/screen_ratio;
set(h_leg,'position', [...
    (marg_l+x_sp*(cm(sp))+gap_w*(cm(sp)-1)-h_pos1(3)), ...
    (marg_b+y_sp*(cm(sp)-1)+gap_h*(cm(sp)-1)), ...
    h_pos1(3), ...
    h_pos1(4)]);

set(ax,'layer','top','color',fig_color,...
    'fontsize',font_size,'tickdir','out',...
    'ticklength',fig_tick_length, ...
    'xtick',115:2:155,'ytick',-20:2:20)
if row_ind(sp) ~= rowN
    set(gca,'xticklabel','')
elseif sp == 4
    xlabel('Longitude')
end
if col_ind(sp) == 1
    ylabel('Transport ($Sv$)')
end

% Cape Leeuwin
hh = arrow([115,-19.7], [115,-17.7], small_arr, ...
    'Facecolor', 'k', 'edgecolor', 'k');
text(115,-19.7,'Cape Leeuwin', 'fontsize', font_size, ...
    'color', 'k','Rotation',-45)
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');
% Albany
hh = arrow([118,-19.7], [118,-17.7], small_arr, ...
    'Facecolor', 'k', 'edgecolor', 'k');
text(118,-19.7,'Albany', 'fontsize', font_size, ...
    'color', 'k','Rotation',-45)
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');
% Cape Pasley
hh = arrow([123.4,-19.7], [123.4,-17.7], small_arr, ...
    'Facecolor', 'k', 'edgecolor', 'k');
text(123.4,-19.7,'Cape Pasley', 'fontsize', font_size, ...
    'color', 'k','Rotation',-45)
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');
% West GAB
hh = arrow([128,-19.7], [128,-17.7], small_arr, ...
    'Facecolor', 'k', 'edgecolor', 'k');
text(128,-19.7,'West GAB', 'fontsize', font_size, ...
    'color', 'k','Rotation',-45)
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');  
% % East GAB
% hh = arrow([132,-19.7], [132,-17.7], small_arr, ...
%     'Facecolor', 'k', 'edgecolor', 'k');
% text(132,-19.7,'East GAB', 'fontsize', font_size, ...
%     'color', 'k','Rotation',-45)
% set(get(get(hh,'Annotation'),'LegendInformation'), ...
%     'IconDisplayStyle','off');  
% Cape Carnot
hh = arrow([135.7,-19.7], [135.7,-17.7], small_arr, ...
    'Facecolor', 'k', 'edgecolor', 'k');
text(135.7,-19.7,'Cape Carnot', 'fontsize', font_size, ...
    'color', 'k','Rotation',-45)
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');  
% Portland
hh = arrow([141.3,-19.7], [141.3,-17.7], small_arr, ...
    'Facecolor', 'k', 'edgecolor', 'k');
text(141.3,-19.7,'Portland', 'fontsize', font_size, ...
    'color', 'k','Rotation',-45)
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');  
%  King Is.
hh = arrow([144,-19.7], [144,-17.7], small_arr, ...
    'Facecolor', 'k', 'edgecolor', 'k');
text(144,-19.7,' King Island', 'fontsize', font_size, ...
    'color', 'k','Rotation',-45)
set(get(get(hh,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle','off');  
% % South East Cape
% hh = arrow([146.7,-19.7], [146.7,-17.7], small_arr, ...
%     'Facecolor', 'k', 'edgecolor', 'k');
% text(146.7,-19.7,'South East C.', 'fontsize', font_size, ...
%     'color', 'k','Rotation',-45)
% set(get(get(hh,'Annotation'),'LegendInformation'), ...
%     'IconDisplayStyle','off'); 


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
close all


%%
play(bird_f_path,[1 (get(bird_f_path, 'SampleRate')*3)]);

