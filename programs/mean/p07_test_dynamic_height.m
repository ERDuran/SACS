%% Inspect maps of dynamic heights
close
figure
set(gcf,'units','normalized','outerposition',[0 0.05 1 0.95])

lat_min = -31;
lat_max = -46;

lon_min = 108;
lon_max = 150;

lat = lat(lat >= lat_max & lat <= lat_min);
lon = lon(lon >= lon_min & lon <= lon_max);

p_ind = pres == p_desired;

dynh_plot = dynh_plot(...
    lat >= lat_max & lat <= lat_min, ...
    lon >= lon_min & lon <= lon_max, ...
    p_ind);

min_dynh_2000 = min(min(dynh_plot));
max_dynh_2000 = max(max(dynh_plot));

levels = min_dynh_2000:increment:max_dynh_2000;

surf(lon,lat,dynh_plot), shading interp

colormap(flipud(othercolor('Spectral9', ...
    length(levels)-1)))

caxis([min_dynh_2000 max_dynh_2000])

aspect_now = pbaspect;
pbaspect(aspect_now.*[1 1 0.2])

axis([lon_min lon_max lat_max lat_min])

colorbar('YTick',levels)


set(gca,'color',[0.7 0.7 0.7])
title(['Dynamic height at ' num2str(pres(p_ind)) ' m derived from CARS T and S'])
grid
font_size = 6;
set(gca,'layer','top', 'tickdir', 'out', ...
    'fontsize',font_size, 'fontweight','bold')


