%% fig 2: meridional cross sections of zonal velocity differences
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_currents'])
load([data_path 'SACS_data/aus8_figures'])
load([data_path 'SACS_data/KDau_rho'])
load([data_path 'SACS_data/aus8_rho'])


%%
close
fig = figure(1);
set(fig,'units','centimeters','paperunits','centimeters', ...
    'inverthardcopy','off',...
    'position',[0 0 20 14]);

subplot(2,1,1)
pcolor(aus8_rho.mean(:,:,1))
caxis([1024.5, 1027])
colormap(othercolor('Spectral10',10))
shading interp
axis ij
colorbar
set(gca, 'color',[0.7 0.7 0.7])
title('Density at the surface in CARS-Aus8')

subplot(2,1,2)
pcolor(KDau_rho.mean(:,:,1))
caxis([1024.5, 1027])
colormap(othercolor('Spectral10',10))
shading interp
axis ij
colorbar
set(gca, 'color',[0.7 0.7 0.7])
title('Density at the surface in MOM01-75z')

print(fig, ...
    ['~/Dropbox/SACS/draft/REVISION/absolute_density'], ...
    '-dpng', '-r300')


fig = figure(2);
set(fig,'units','centimeters','paperunits','centimeters', ...
    'inverthardcopy','off',...
    'position',[0 0 20 7]);
pcolor(aus8_rho.mean(:,:,1) - KDau_rho.mean(:,:,1))
caxis([-0.5, 0.5])
colormap(othercolor('RdBu11',10))
shading interp
axis ij
colorbar
set(gca, 'color',[0.7 0.7 0.7])
title('Density difference at the surface: CARS-Aus8 - MOM01-75z')

print(fig, ...
    ['~/Dropbox/SACS/draft/REVISION/density_anomaly'], ...
    '-dpng', '-r300')