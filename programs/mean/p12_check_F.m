%% test divergence = 0
% clc
path(pathdef)

% set up main directory
cd /home/z5100174/Desktop/MATLAB

% Add path to the data to be loaded
addpath(genpath('functions'))

% Add path to the data to be loaded
addpath cars_out

% Add path to the gsw TEOS-10 library (including subfolders)
addpath(genpath('teos10_library'))

clear 
% vn is a function that can store the name of a variable
% hence vn(x) returns 'x' as a string.
vn = @(x) inputname(1);

% load geostrophic velocity
load aus8_ZD_method

lat = aus8_ZD_method.lat_F;
lon = aus8_ZD_method.lon_F;
F = aus8_ZD_method.F;


%% plot map one order
close
figure
x_pos = 0; % x figure position
y_pos = 0.05; % y figure position
x_lgth = 0.95; % x figure length
y_lgth = 0.99; % y figure length
set(gcf,'units','normalized','outerposition', ...
    [x_pos y_pos x_lgth y_lgth]) % apply onto current figure

rowN = 2; % number of rows
colN = 1; % number of columns
gap_w = 0.04; % gap width between subplots
gap_h = 0.08; % gap height between subplots
marg_b = 0.05; % bottom margin
marg_t = 0.04; % top margin
marg_l = 0.05; % left margin
marg_r = 0.05; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length
h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]); % creates subplots

font_size = 8; % font size
background_c = [0.7 0.7 0.7]; % background colour 

order = 6; % cbar order of magnitude

axes(h_axes_sp(1)) % call one subplot

pcolor(lon,lat,F), shading interp
caxis([-10^-order 10^-order])
colormap(othercolor('RdGy11'))
colorbar
set(gca,'color',background_c)
title(sprintf('F = div(U,V)'))
grid
set(gca,'layer','top', 'tickdir', 'out', ...
    'fontsize',font_size, 'fontweight','bold')


order = 17; % cbar order of magnitude

axes(h_axes_sp(2)) % call one subplot

pcolor(lon,lat,F), shading interp
caxis([-10^-order 10^-order])
colormap(othercolor('RdGy11'))
colorbar
set(gca,'color',background_c)
title(sprintf('F = div(U,V)'))
grid
set(gca,'layer','top', 'tickdir', 'out', ...
    'fontsize',font_size, 'fontweight','bold')


export_fig('figures/p3_F_tests/F', '-m2', '-nocrop');

