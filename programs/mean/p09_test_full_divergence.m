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
load aus8_geostrophy

% load topography mask
load topog_mask

lat = aus8_geostrophy.div.lat;
lon = aus8_geostrophy.div.lon;
pres = aus8_geostrophy.pres;

div_0 = aus8_geostrophy.div.div_0;
div_2000_HH = aus8_geostrophy.div.div_2000_HH;


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
colN = 3; % number of columns
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

for jj = 1 : colN % go through columns
    
    % set pressure level
    if jj == 1
        pres_now = 0;
    elseif jj == 2
        pres_now = 1000;
    elseif jj == 3
        pres_now = 2000;
    end
    pres_ind = pres == pres_now;
    
    axes(h_axes_sp(jj)) % call one subplot
    
    pcolor(lon,lat,div_0(:,:,pres_ind)), shading interp
    caxis([-10^-order 10^-order])
    colormap(othercolor('RdGy11'))
    colorbar
    set(gca,'color',background_c)
    title(sprintf(['Horiz. Divergence at ' num2str(pres_now) ...
        ' dbars from Geostrophic Velocities \n' ...
        'derived on an f-plane ' ...
        'wrt the surface']))
    grid
    set(gca,'layer','top', 'tickdir', 'out', ...
        'fontsize',font_size, 'fontweight','bold')
    
    axes(h_axes_sp(jj+3)) % call one subplot
    
    pcolor(lon,lat,div_2000_HH(:,:,pres_ind)), shading interp
    caxis([-10^-order 10^-order])
    colormap(othercolor('RdGy11'))
    colorbar
    set(gca,'color',background_c)
    title(sprintf(['Horiz. Divergence  at ' num2str(pres_now) ...
        ' dbars from Geostrophic Velocities \n' ...
        'derived on an f-plane ' ...
        'wrt 2000 bars or the sea bottom']))
    grid
    set(gca,'layer','top', 'tickdir', 'out', ...
        'fontsize',font_size, 'fontweight','bold')
    
    
end

export_fig(['figures/div_tests/order_' num2str(order)], '-m2', '-nocrop');


% plot map the other order
close
figure
x_pos = 0; % x figure position
y_pos = 0.05; % y figure position
x_lgth = 0.95; % x figure length
y_lgth = 0.99; % y figure length
set(gcf,'units','normalized','outerposition', ...
    [x_pos y_pos x_lgth y_lgth]) % apply onto current figure

rowN = 2; % number of rows
colN = 3; % number of columns
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

order = 20; % cbar order of magnitude

for jj = 1 : colN % go through columns
    
    % set pressure level
    if jj == 1
        pres_now = 0;
    elseif jj == 2
        pres_now = 1000;
    elseif jj == 3
        pres_now = 2000;
    end
    pres_ind = pres == pres_now;
    
    axes(h_axes_sp(jj)) % call one subplot
    
    pcolor(lon,lat,div_0(:,:,pres_ind)), shading interp
    caxis([-10^-order 10^-order])
    colormap(othercolor('RdGy11'))
    colorbar
    set(gca,'color',background_c)
    title(sprintf(['Horiz. Divergence at ' ...
        num2str(pres_now) ' dbars from Geostrophic Velocities \n' ...
        'derived on an f-plane ' ...
        'wrt the surface']))
    grid
    set(gca,'layer','top', 'tickdir', 'out', ...
        'fontsize',font_size, 'fontweight','bold')
    
    axes(h_axes_sp(jj+3)) % call one subplot
    
    pcolor(lon,lat,div_2000_HH(:,:,pres_ind)), shading interp
    caxis([-10^-order 10^-order])
    colormap(othercolor('RdGy11'))
    colorbar
    set(gca,'color',background_c)
    title(sprintf(['Horiz. Divergence at ' ...
        num2str(pres_now) ' dbars from Geostrophic Velocities \n' ...
        'derived on an f-plane ' ...
        'wrt 2000 bars or the sea bottom']))
    grid
    set(gca,'layer','top', 'tickdir', 'out', ...
        'fontsize',font_size, 'fontweight','bold')
    
    
end

export_fig(['figures/p1_div_tests/order_' num2str(order)], '-m2', '-nocrop');
