%% load data
% note since I want to plot a map using quiver,
% u and v must be on the same grid (T and S grid)
% therefore this will ruin the zero divergence property
% which is ok here because we just want to look at maps of
% geostrophic velocity fields


[u_interp2, v_interp2] = ...
    deal(NaN(length(lat_u),length(lon_v),length(pres_mid)));
for kk = 1 : length(pres_mid)
    u_interp2(:,:,kk) = interp2(lon_u,lat_u,u_plot(:,:,kk),lon_v,lat_u);
    v_interp2(:,:,kk) = interp2(lon_v,lat_v,v_plot(:,:,kk),lon_v,lat_u);
end

speed = sqrt(u_interp2.^2 + v_interp2.^2);


%%
% Xsects
merid_Xsect_lat = [ ...
    -34; %1
    -35; %2
    -35; %3
    -35; %4
    -34; %5
    -34; %6
    -34; %7
    -34; %8
    -34; %9
    -34; %10
    -33; %11
    -33; %12
    -32; %13
    -32; %14
    -32; %15
    -32; %16
    -32; %17
    -32; %18
    -33; %19
    -33; %20
    -34; %21
    -35; %22
    -36; %23
    -36; %24
    -37; %25
    -37; %26
    -38; %27
    -38; %28
    -39; %29
    -40; %30
    -42; %31
    -43; %32
    -43; %33
    -43; %35
    ]';

merid_Xsect_lat(2,:) = merid_Xsect_lat(1,:) - 5;


merid_Xsect_lon = [115:148; 115:148];

fig1 = figure(1);

set(gcf,'units','normalized','outerposition',[0 0.1 0.99 0.97])

rowN = 4; colN = 9;
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (flipud(repmat((1:rowN)',1,colN)))';
gap_w = 0.015; % gap width between subplots
gap_h = 0.04; % gap height between subplots
marg_b = 0.03; % bottom margin
marg_t = 0.02; % top margin
marg_l = 0.02; % left margin
marg_r = 0.005; % right margin
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
    
    
    lat_ind = lat_u <= merid_Xsect_lat(1,ss) & lat_u >= merid_Xsect_lat(2,ss);
    
    lon_ind = lon_u <= merid_Xsect_lon(1,ss) & lon_u >= merid_Xsect_lon(2,ss);
    
    lat_now = lat_u(lat_ind);
    
    [CMAP,LEV,WID,CLIM] = cmjoin({...
        flipud(othercolor('Blues7')),othercolor('Reds7')}, ...
        [-0.1 0 0.1], [0.01 0.01]);
    cmap = colormap(h_axes_sp(1), CMAP);
    
    pcol = pcolor(lat_now, pres_mid', squeeze(u_plot(lat_ind,lon_ind,:))');
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
    
    set(gca,'ytick',0:200:2000,'yticklabel','')
    if ss == 1 || ss == 10 || ss == 19 || ss == 28
        set(gca, 'yticklabel',0:200:2000)
    end
    
    set(gca,'xtick',merid_Xsect_lat(2,ss):merid_Xsect_lat(1,ss))
    
    grid
    font_size = 8;
    set(gca,'layer','top','color',[0.7 0.7 0.7],'fontsize',font_size)
    
end


%%
% fig2 = figure(2);
% 
% set(gcf,'units','normalized','outerposition',[0 0.1 0.99 0.97])
% 
% rowN = 4; colN = 9;
% col_ind = (repmat(1:colN,rowN,1))';
% row_ind = (flipud(repmat((1:rowN)',1,colN)))';
% gap_w = 0.015; % gap width between subplots
% gap_h = 0.04; % gap height between subplots
% marg_b = 0.03; % bottom margin
% marg_t = 0.02; % top margin
% marg_l = 0.02; % left margin
% marg_r = 0.005; % right margin
% x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
% y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length
% 
% h_axes_sp = tight_subplot(rowN, colN, ...
%     [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);
% 
% for ss = 1 : length(merid_Xsect_lat) + 1
%     axes(h_axes_sp(ss))
%     
%     if ss == length(merid_Xsect_lat) + 1
%         [CMAP,LEV,WID,CLIM] = cmjoin({...
%             flipud(othercolor('Blues7')),othercolor('Reds7')}, ...
%             [-0.1 0 0.1], [0.01 0.01]);
%         colorbar
%         axis off
%         break
%     end
%     
%     
%     lat_ind = lat_v <= merid_Xsect_lat(1,ss) & lat_v >= merid_Xsect_lat(2,ss);
%     
%     lon_ind = abs(lon_v-merid_Xsect_lon(1,ss)) == ...
%         min(abs(lon_v-merid_Xsect_lon(1,ss)));
%     
%     lon_now = lon_v(lon_ind);
%     lat_now = lat_v(lat_ind);
%     
%     [CMAP,LEV,WID,CLIM] = cmjoin({...
%         flipud(othercolor('Blues7')),othercolor('Reds7')}, ...
%         [-0.1 0 0.1], [0.01 0.01]);
%     cmap = colormap(h_axes_sp(1), CMAP);
%     
%     pcol = pcolor(lat_now, pres_mid', squeeze(v_plot(lat_ind,lon_ind,:))');
%     axis ij
%     shading interp
% 
%     caxis([-0.1 0.1])
%     
%     hold on
%     
%     
%     % h_quiv = quivers(...
%     %     lon_mg(1:nn:end,1:nn:end), lat_mg(1:nn:end,1:nn:end), ...
%     %     u_interp2(1:nn:end,1:nn:end), v_interp2(1:nn:end,1:nn:end), ...
%     %     quiv_S, 1, 'm/s', 'k');
%     
%     
%     title(['u (m/s) at ' num2str(merid_Xsect_lon(1,ss))...
%         ' ^\circE.'])
%     
%     set(gca,'ytick',0:200:2000,'yticklabel','')
%     if ss == 1 || ss == 10 || ss == 19 || ss == 28
%         set(gca, 'yticklabel',0:200:2000)
%     end
%     
%     set(gca,'xtick',merid_Xsect_lat(2,ss):merid_Xsect_lat(1,ss))
%     
%     grid
%     font_size = 8;
%     set(gca,'layer','top','color',[0.7 0.7 0.7],'fontsize',font_size)
%     
% end


%%
fig3 = figure(3);

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


pcolor(lon_v, lat_u, speed_s);
shading interp
colorbar

caxis([0 0.04])

hold on

nn = 2;
quiv_S = 4;

[lon_mg, lat_mg]=meshgrid(lon_v,lat_u);

quivers(...
    lon_mg(1:nn:end,1:nn:end), lat_mg(1:nn:end,1:nn:end), ...
    u_interp2_s(1:nn:end,1:nn:end), v_interp2_s(1:nn:end,1:nn:end), ...
    quiv_S, 1, 'm/s', 'k');

title(['Geostrophic velocity (m/s) field at 5 db.' ...
    ' Location of the meridional sections'])

hold on
plot(merid_Xsect_lon,merid_Xsect_lat,'r','linewidth',2)

grid 
font_size = 8;
set(gca,'layer','top','color',[0.7 0.7 0.7],'fontsize',font_size)

axes(h_axes_sp(2))

colormap(h_axes_sp(2), flipud(othercolor('Blues9')));

d_ind = pres_mid == pres_mid_2;

pcolor(lon_v, lat_u, speed(:,:,d_ind));
shading interp
colorbar

caxis([0 0.04])

hold on

[lon_mg, lat_mg]=meshgrid(lon_v,lat_u);

h_quiv = quivers(...
    lon_mg(1:nn:end,1:nn:end), lat_mg(1:nn:end,1:nn:end), ...
    u_interp2(1:nn:end,1:nn:end,d_ind), v_interp2(1:nn:end,1:nn:end,d_ind), ...
    quiv_S, 1, 'm/s', 'k');

title(['Geostrophic velocity (m/s) field at 405 db.' ...
    ' Location of the meridional sections'])

hold on
plot(merid_Xsect_lon,merid_Xsect_lat,'r','linewidth',2)

grid 
font_size = 8;
set(gca,'layer','top','color',[0.7 0.7 0.7],'fontsize',font_size)

