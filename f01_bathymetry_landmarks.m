%% fig 2: map bathymetry and landmarks. In black and white
clearvars('-except', '*_path')

load topog_mask

lat = topog_mask.etopo5.lat;
lon = topog_mask.etopo5.lon;
topog_smoothed_raw = topog_mask.etopo5.topog_smoothed;


%% plot map
% 1) figure set-up
close all
font_size = 6;
fig1 = figure(1);
% magn = 3;
% set(gcf,'PaperType','A4', ...
%     'paperOrientation', 'portrait', ...
%     'paperunits','centimeters', ...
%     'PaperPosition',[0, 0, (21-4)/2, (29.7-4)]*magn);
set(gcf, 'units', 'normalized', 'outerposition', [0 0 0.4 0.4]);
rowN = 1; colN = 1;
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = - 0.04; % gap width between subplots
gap_h = 0.04; % gap height between subplots
marg_b = 0.05; % bottom margin
marg_t = 0.05; % top margin
marg_l = 0.04; % left margin
marg_r = 0.02; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length
h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

% 2) all data set-up
sp = 1;
axes(h_axes_sp(sp))
topog_smoothed = topog_smoothed_raw;
topog_smoothed(isnan(topog_smoothed_raw)) = 0;

% 3) contourf set-up
contourf_contours = [0 0];
cmap_levels = length(contourf_contours);
cmap = [0.8 0.8 0.8];
colormap(h_axes_sp(sp), cmap);

% 4) plot contourf
contourf(...
    lon, ...
    lat, ...
    topog_smoothed, ...
    contourf_contours);
hold on
caxis([contourf_contours(1) contourf_contours(end)]);
axis([112 150 -45 -31])


% 5) Theta contours set-up
contour_contours = -[100 200 500 1000 1500 2000];

% 6) plot Theta labelled contours
[c, h] = contour(...
    lon, ...
    lat, ...
    topog_smoothed, ...
    contour_contours, ...
    'k', ...
    'linewidth',0.1);
clabel(c,h, 'manual',...
    'fontsize', font_size)


% 7) prepare text
% 8) place landmarks
% Cape Leeuwin
text(115,-34.15,'Cape Leeuwin', 'fontsize', font_size)
% Cape Pasley
text(123.4,-33.85,'Cape Pasley', 'fontsize', font_size)
% Great Australian Bight
text(127.5,-32.65,'Great Australian Bight', 'fontsize', font_size)
% Cape Carnot
text(135.5,-34.85,'Cape Carnot', 'fontsize', font_size)
% Spencer Gulf
text(136.9,-33.8,'Spencer Gulf', 'fontsize', font_size)
% St. Vincent Gulf
text(138.1,-34.24,'St. Vincent Gulf', 'fontsize', font_size)
% Portland
text(141.3,-38.33,'Portland', 'fontsize', font_size)
% Cape Otway
text(143.5,-38.75,'Cape Otway', 'fontsize', font_size)
% South East Cape
text(146.7,-43.62,'South East Cape', 'fontsize', font_size)

% 9) title, grid, background and fonts
title('ETOPO5 bathymetry and landmarks')
grid
set(gca,'layer','top','color',[1 1 1],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end


% Save
% outputls = ls(figures_path);
% scriptname = mfilename;
% if ~contains(outputls, scriptname)
%     mkdir(figures_path, scriptname)
% end
% export_fig(fig1, [figures_path mfilename '/' scriptname(1:3) '_'], ...
%     '-m3')
% close

