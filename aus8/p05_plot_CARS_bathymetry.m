%% Inspect maps of dynamic heights
clc
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

load aus8_coast_geostrophy
load aus8_coast_TEOS10
load topog_mask
load Helland_Hansen


%% q
lat = aus8_coast_TEOS10.lat;
lon = aus8_coast_TEOS10.lon;
pres = aus8_coast_TEOS10.pres;

sigma = aus8_coast_TEOS10.sigma.mean;

CARS_bathymetry = ones(length(lat), length(lon)) * 10;
for ii = 1 : length(lat)
    for jj = 1 : length(lon)
        bottom_ind = find(isfinite(sigma(ii,jj,:)),1,'last');
        
        if isempty(bottom_ind)
            continue
        else
            CARS_bathymetry(ii,jj) = -pres(bottom_ind);
        end
    end
end

etopo5 = topog_mask.etopo5.topog_smoothed;
etopo5(etopo5<-2000) = -2000;
etopo5(isnan(etopo5)) = 10;
lon_etopo5 = topog_mask.etopo5.lon;
lat_etopo5 = topog_mask.etopo5.lat;
smoothing_window = topog_mask.etopo5.smoothing_window;

topog = topog_mask.topog.topog;
lat_topog = topog_mask.topog.lat;
lon_topog = topog_mask.topog.lon;


%%
close
figure
set(gcf,'units','normalized','outerposition',[0 0.05 1 0.95])

rowN = 2; colN = 2;
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (flipud(repmat((1:rowN)',1,colN)))';
gap_w = 0.01; % gap width between subplots
gap_h = 0.06; % gap height between subplots
marg_b = 0.05; % bottom margin
marg_t = 0.04; % top margin
marg_l = 0.01; % left margin
marg_r = 0.01; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length

h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

axes(h_axes_sp(1))

lat_min = -31;
lat_max = -46;

lon_min = 112;
lon_max = 152;

lat = aus8_coast_geostrophy.lat;
lon = aus8_coast_geostrophy.lon;
bathymetry = CARS_bathymetry;

min_CARS_bathymetry = -2000;
max_CARS_bathymetry = 0;

surf(lon,lat,bathymetry), %shading flat

z_magn = 0.6;

[CMAP,LEV,WID,CLIM] = cmjoin({flipud(...
    othercolor('RdGy9')), [0.4039,0,0.1215], othercolor('Spectral11'), [0 0 0]}, ...
    [-2050 -250 -195 5 15], ...
    [100 55 10 10]);
colormap(CMAP)

axis([lon_min lon_max lat_max lat_min])

aspect_now = pbaspect;
pbaspect(aspect_now.*[1 1 z_magn])

set(gca,'color',[0 0 0])
title('Original CARS bathymetry map (meters). Rainbow shading increment: 10 m')
grid on
set(gca,'GridColor',[1 1 1], 'GridAlpha', 0.4, 'GridLineStyle', '--');
font_size = 8;
set(gca,'layer','top', 'tickdir', 'out', ...
    'fontsize',font_size, 'fontweight','bold')


axes(h_axes_sp(2))

lat = lat_etopo5;
lon = lon_etopo5;
bathymetry = etopo5;


surf(lon,lat,bathymetry), shading flat

colormap(CMAP)

axis([lon_min lon_max lat_max lat_min])
aspect_now = pbaspect;
pbaspect(aspect_now.*[1 1 z_magn])

set(gca,'color',[0 0 0])
title(['ETOPO-5 bathymetry map (meters) smoothed out within a window of ±' ...
    num2str(smoothing_window*0.025) '\circ lat and lon'])
grid on
set(gca,'GridColor',[1 1 1], 'GridAlpha', 0.4, 'GridLineStyle', '--');
font_size = 8;
set(gca,'layer','top', 'tickdir', 'out', ...
    'fontsize',font_size, 'fontweight','bold')

axes(h_axes_sp(3))

lat = lat_topog;
lon = lon_topog;
bathymetry = topog;

surf(lon,lat,bathymetry), %shading flat

colormap(CMAP)

axis([lon_min lon_max lat_max lat_min])
aspect_now = pbaspect;
pbaspect(aspect_now.*[1 1 z_magn])

set(gca,'color',[0 0 0])
title(['Smoothed ETOPO-5 map interpolated onto CARS grid (meters).', ...
    ' Linear interpolation used.'])
grid on
set(gca,'GridColor',[1 1 1], 'GridAlpha', 0.4, 'GridLineStyle', '--');
font_size = 8;
set(gca,'layer','top', 'tickdir', 'out', ...
    'fontsize',font_size, 'fontweight','bold')

axes(h_axes_sp(4))

h_surf = surf(lon,lat,bathymetry); shading flat

colormap(CMAP)

axis([lon_min lon_max lat_max lat_min])
aspect_now = pbaspect;
pbaspect(aspect_now.*[1 1 z_magn])

set(gca,'color',[0 0 0])
title(['Smoothed ETOPO-5 map interpolated onto CARS grid (meters).', ...
    ' Area zoomed in.'])
grid on
set(gca,'GridColor',[1 1 1], 'GridAlpha', 0.4, 'GridLineStyle', '--');
font_size = 8;
set(gca,'layer','top', 'tickdir', 'out', ...
    'fontsize',font_size, 'fontweight','bold')

colorbar('horizontal','ytick',[-2000:100:0])
caxis(CLIM)


%%
ax_list = {...
    [lon_min lon_max lat_max lat_min -2000 0], ...
    [113 120 -36 -31 -2000 10], ...
    [120 132 -36 -31 -2000 0], ...
    [132 143 -41 -31 -2000 0], ...
    [141 146 -45 -36 -2000 0], ...
    [146 150 -46 -37 -2000 0]};

for pp = 2%1 : 6
    if pp == 1
        pic = 'p1_all';
        ax = ax_list{pp};
        v1 = 0; v2 = 90;
        
    elseif pp == 2
        pic = 'p2_southwest';
        ax = ax_list{pp};
        v1 = -5; v2 = 35;
        
    elseif pp == 3
        pic = 'p3_GAB';
        ax = ax_list{pp};
        v1 = -20; v2 = 25;
        
    elseif pp == 4
        pic = 'p4_ADL';
        ax = ax_list{pp};
        v1 = -32; v2 = 30;
        
    elseif pp == 5
        pic = 'p5_west_TAS';
        ax = ax_list{pp};
        v1 = -5; v2 = 60;
        
    elseif pp == 6
        pic = 'p6_east_TAS';
        ax = ax_list{pp};
        v1 = 40; v2 = 20;
    end
    
    for sp = 1 : 3
        axes(h_axes_sp(sp))
        if sp == 1 || sp == 3
            if strcmp(pic,'p1_all')
                shading flat
            else
                shading faceted
            end
        end
        
        axis(ax)
        view_now = view(v1,v2);
    end
    
    if ~strcmp(pic,'p1_all')
        x_axis = [ax_list{pp}(1), ax_list{pp}(2), ...
            ax_list{pp}(2), ax_list{pp}(1), ax_list{pp}(1)];
        y_axis = [ax_list{pp}(3), ax_list{pp}(3), ...
            ax_list{pp}(4), ax_list{pp}(4), ax_list{pp}(3)];
        
        axes(h_axes_sp(4))
        
        if exist('h_plot')
            set(h_plot,'Visible','off')
        end
        
        hold on
        h_plot = plot(x_axis, y_axis, 'b', 'linewidth', 3);
    end
    
    
    export_fig(['figures/bathymetry/smooth' ...
        num2str(smoothing_window) '/' pic], '-m3', '-nocrop');
end


%%
% save('cars_out/views', ...
%     'views')
% disp(['views saved in ' ...
%     'cars_out/views !'])
