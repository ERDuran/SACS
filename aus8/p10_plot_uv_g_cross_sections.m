%% load data
% note since I want to plot a map using quiver,
% u and v must be on the same grid (T and S grid)
% therefore this will ruin the zero divergence property
% which is ok here because we just want to look at maps of
% geostrophic velocity fields


[u_interp2, v_interp2] = deal(NaN(length(lat_v),length(lon_u),length(pres_mid)));
for kk = 1 : length(pres_mid)
    u_interp2(:,:,kk) = interp2(lon_u,lat_u,u_g_mid(:,:,kk),lon_u,lat_v);
    v_interp2(:,:,kk) = interp2(lon_v,lat_v,v_g_mid(:,:,kk),lon_u,lat_v);
end

speed = sqrt(u_interp2.^2 + v_interp2.^2);


%%
% Xsects
merid_Xsect_lat = [ ...
    -34 -39; %1
    -34 -39; %2
    -34 -39; %3
    -34 -39; %4
    -34 -39; %5
    -33 -38; %6
    -32 -37; %7
    -31 -36; %8
    -32 -37; %9
    -33 -38; %10
    -34 -39; %11
    -36 -41; %12
    -36 -41; %13
    -38 -43; %14
    -38 -43]';
merid_Xsect_lon = [ ...
    115 115; %1
    117 117; %2
    119 119; %3
    121 121; %4
    123 123; %5
    125 125; %6
    127 127; %7
    129 129; %8
    131 131; %9
    133 133; %10
    135 135; %11
    137 137; %12
    139 139; %13
    141 141; %14
    143 143]';

fig1 = figure(1);

set(gcf,'units','normalized','outerposition',[0 0.1 0.9 0.9])

rowN = 2; colN = 8;
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (flipud(repmat((1:rowN)',1,colN)))';
gap_w = 0.02; % gap width between subplots
gap_h = 0.06; % gap height between subplots
marg_b = 0.05; % bottom margin
marg_t = 0.02; % top margin
marg_l = 0.05; % left margin
marg_r = 0.05; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length

h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

for ss = 1 : length(merid_Xsect_lat) + 1
    axes(h_axes_sp(ss))
    
    if ss == length(merid_Xsect_lat) + 1
        [CMAP,LEV,WID,CLIM] = cmjoin({...
            flipud(othercolor('Blues7')),othercolor('Reds7')}, ...
            [-0.1 0 0.1], [0.01 0.01]);
        colorbar
        axis off
        break
    end
    
    
    lat_ind = lat_v <= merid_Xsect_lat(1,ss) & lat_v >= merid_Xsect_lat(2,ss);
    
    lon_ind = lon_u <= merid_Xsect_lon(1,ss) & lon_u >= merid_Xsect_lon(2,ss);
    
    lat_now = lat_v(lat_ind);
    lon_now = lon_u(lon_ind);
    
    [CMAP,LEV,WID,CLIM] = cmjoin({...
        flipud(othercolor('Blues7')),othercolor('Reds7')}, ...
        [-0.1 0 0.1], [0.01 0.01]);
    cmap = colormap(h_axes_sp(1), CMAP);
    
    pcol = pcolor(lat_now, pres_mid', squeeze(u_g_mid(lat_ind,lon_ind,:))');
    axis ij
    shading interp

    caxis([-0.1 0.1])
    
    hold on
    
    
    % h_quiv = quivers(...
    %     lon_mg(1:nn:end,1:nn:end), lat_mg(1:nn:end,1:nn:end), ...
    %     u_interp2(1:nn:end,1:nn:end), v_interp2(1:nn:end,1:nn:end), ...
    %     quiv_S, 1, 'm/s', 'k');
    
    
    title(['u (m/s) at ' num2str(merid_Xsect_lon(1,ss))...
        ' ^\circE.'])
    
    set(gca,'yticklabel','')
    if ss == 1 || ss == 5
        set(gca,'yticklabel',0:200:2000)
    end
      
    grid
    font_size = 8;
    set(gca,'layer','top','color',[0.7 0.7 0.7],'fontsize',font_size)
    
end






%%
fig2 = figure(2);

set(gcf,'units','normalized','outerposition',[0 0.1 0.6 0.99])

rowN = 2; colN = 1;
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (flipud(repmat((1:rowN)',1,colN)))';
gap_w = 0.02; % gap width between subplots
gap_h = 0.06; % gap height between subplots
marg_b = 0.05; % bottom margin
marg_t = 0.02; % top margin
marg_l = 0.05; % left margin
marg_r = 0.05; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length

h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

axes(h_axes_sp(1))

colormap(h_axes_sp(1), flipud(othercolor('Blues9')));

d_ind = pres_mid == pres_mid_1;

u_interp2_s = u_interp2(:,:,d_ind);
v_interp2_s = v_interp2(:,:,d_ind);
speed_s = speed(:,:,d_ind);

for ii = 1 : 10
    max_speed = find(speed == max(max(speed_s)));
    
    [u_interp2_s(max_speed), v_interp2_s(max_speed), speed_s(max_speed)] = ...
        deal(NaN);
end


pcolor(lon_u, lat_v, speed_s);
shading interp
colorbar

caxis([0 0.04])

hold on

nn = 2;
quiv_S = 4;

[lon_mg, lat_mg]=meshgrid(lon_u,lat_v);

quivers(...
    lon_mg(1:nn:end,1:nn:end), lat_mg(1:nn:end,1:nn:end), ...
    u_interp2_s(1:nn:end,1:nn:end), v_interp2_s(1:nn:end,1:nn:end), ...
    quiv_S, 1, 'm/s', 'k');

title(['Geostrophic velocity (m/s) field at the surface.' ...
    ' Location of the meridional sections'])

hold on
plot(merid_Xsect_lon,merid_Xsect_lat,'r','linewidth',2)

grid 
font_size = 8;
set(gca,'layer','top','color',[0.7 0.7 0.7],'fontsize',font_size)

axes(h_axes_sp(2))

colormap(h_axes_sp(2), flipud(othercolor('Blues9')));

d_ind = pres_mid == pres_mid_2;

pcolor(lon_u, lat_v, speed(:,:,d_ind));
shading interp
colorbar

caxis([0 0.04])

hold on

[lon_mg, lat_mg]=meshgrid(lon_u,lat_v);

h_quiv = quivers(...
    lon_mg(1:nn:end,1:nn:end), lat_mg(1:nn:end,1:nn:end), ...
    u_interp2(1:nn:end,1:nn:end,d_ind), v_interp2(1:nn:end,1:nn:end,d_ind), ...
    quiv_S, 1, 'm/s', 'k');

title(['Geostrophic velocity (m/s) field at 400 db.' ...
    ' Location of the meridional sections'])

hold on
plot(merid_Xsect_lon,merid_Xsect_lat,'r','linewidth',2)

grid 
font_size = 8;
set(gca,'layer','top','color',[0.7 0.7 0.7],'fontsize',font_size)

