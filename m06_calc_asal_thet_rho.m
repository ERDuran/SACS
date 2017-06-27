%% Calculate conservative temperature Theta, absolute salinity asal,
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_temp'])
load([data_path 'SACS_data/aus8_salt'])
load([data_path 'SACS_data/ETOPO5'])

lat = aus8_coor.lat;
lon = aus8_coor.lon;
Months = aus8_coor.Months;
depth = aus8_coor.depth;
pres = gsw_p_from_z(depth,-40);
aus8_coor.pres = pres;
save([data_path 'SACS_data/aus8_coor'], 'aus8_coor')


%% Calculate conservative temperature, absolute salinity
% make templates
[aus8_asal.mean, aus8_thet.mean, aus8_rho.mean] = ...
    deal(NaN([length(lat) length(lon) length(depth)]));
for t = 1 : 12
    aus8_asal.(Months{t}) = ...
        NaN([length(lat), length(lon), length(depth)]);
    aus8_thet.(Months{t}) = ...
        NaN([length(lat), length(lon), length(depth)]);
    aus8_rho.(Months{t}) = ...
        NaN([length(lat), length(lon), length(depth)]);
end

% Calculate using functions from Gibbs seawater toolbox TEOS-10
for p = 1 : length(pres)
    % Mean absolute salinity and conservative temperature from
    % practical salinity and in-situ temperature
    aus8_asal.mean(:,:,p) = gsw_SA_from_SP(...
        aus8_salt.mean(:,:,p),pres(p),lon,lat);
    aus8_thet.mean(:,:,p) = gsw_CT_from_t(...
        aus8_asal.mean(:,:,p),...
        aus8_temp.mean(:,:,p),pres(p));
    
    for t = 1 : 12
        aus8_asal.(Months{t})(:,:,p) = gsw_SA_from_SP(...
            aus8_salt.(Months{t})(:,:,p),pres(p),lon,lat);
        
        aus8_thet.(Months{t})(:,:,p) = gsw_CT_from_t(...
            aus8_asal.(Months{t})(:,:,p),...
            aus8_temp.(Months{t})(:,:,p),pres(p));
    end
    fprintf('p = %4.0f \n',depth(p))
end


% Calculate density
for p = 1 : length(pres)
    aus8_rho.mean(:,:,p) = gsw_rho(...
        aus8_asal.mean(:,:,p),...
        aus8_thet.mean(:,:,p),pres(p));
    
    for t = 1 : 12
        aus8_rho.(Months{t})(:,:,p) = gsw_rho(...
            aus8_asal.(Months{t})(:,:,p),...
            aus8_thet.(Months{t})(:,:,p),pres(p));
    end
end


%% figure set-up mean
close all
font_size = 10;
fig_n = 1;
fig = figure(fig_n);
rowN = 3; colN = 3;
set(gcf,'units','centimeters','position',[0 0 colN*10 rowN*7]);

col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = 0.02; % gap width between subplots
gap_h = 0.04; % gap height between subplots
marg_b = 0.04; % bottom margin
marg_t = 0.04; % top margin
marg_l = 0.04; % left margin
marg_r = 0.02; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length
h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

% 1)
sp = 1;
axes(h_axes_sp(sp))
depth_lvl = 0;
depth_lvl_ind = depth == depth_lvl;
data = aus8_thet.mean(:,:,depth_lvl_ind);
data_min = 7;
data_max = 21;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 mean $\Theta$ $z=' num2str(depth_lvl) '$ ($^{\circ}C$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 2)
sp = 2;
axes(h_axes_sp(sp))
depth_lvl = 0;
depth_lvl_ind = depth == depth_lvl;
data = aus8_asal.mean(:,:,depth_lvl_ind);
data_min = 34;
data_max = 37;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('PuOr8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 mean $S_{A}$ $z=' num2str(depth_lvl) '$ ($g/kg$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 3)
sp = 3;
axes(h_axes_sp(sp))
depth_lvl = 0;
depth_lvl_ind = depth == depth_lvl;
data = aus8_rho.mean(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = 1026.8;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('BrBG8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 mean $\rho$ $z=' num2str(depth_lvl) '$ ($kg/m^3$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 4)
sp = 4;
axes(h_axes_sp(sp))
depth_lvl = -400;
depth_lvl_ind = depth == depth_lvl;
data = aus8_thet.mean(:,:,depth_lvl_ind);
data_min = 7;
data_max = 11;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 mean $\Theta$ $z=' num2str(depth_lvl) '$ ($^{\circ}C$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 5)
sp = 5;
axes(h_axes_sp(sp))
depth_lvl = -400;
depth_lvl_ind = depth == depth_lvl;
data = aus8_asal.mean(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = 35.2;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('PuOr8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 mean $S_{A}$ $z=' num2str(depth_lvl) '$ ($g/kg$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 6)
sp = 6;
axes(h_axes_sp(sp))
depth_lvl = -400;
depth_lvl_ind = depth == depth_lvl;
data = aus8_rho.mean(:,:,depth_lvl_ind);
data_min = 1028.5;
data_max = 1028.8;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('BrBG8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 mean $\rho$ $z=' num2str(depth_lvl) '$ ($kg/m^3$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 7)
sp = 7;
axes(h_axes_sp(sp))
depth_lvl = -1000;
depth_lvl_ind = depth == depth_lvl;
data = aus8_thet.mean(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 mean $\Theta$ $z=' num2str(depth_lvl) '$ ($^{\circ}C$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 8)
sp = 8;
axes(h_axes_sp(sp))
depth_lvl = -1000;
depth_lvl_ind = depth == depth_lvl;
data = aus8_asal.mean(:,:,depth_lvl_ind);
data_min = 34.5;
data_max = 34.62;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('PuOr8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 mean $S_{A}$ $z=' num2str(depth_lvl) '$ ($g/kg$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 9)
sp = 9;
axes(h_axes_sp(sp))
depth_lvl = -1000;
depth_lvl_ind = depth == depth_lvl;
data = aus8_rho.mean(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('BrBG8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 mean $\rho$ $z=' num2str(depth_lvl) '$ ($kg/m^3$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% Save
outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
export_fig(fig, ...
    [figures_path mfilename '/' scriptname(1:3) ...
    '_fig' num2str(fig_n) '_'], ...
    '-m3')
close


%% figure set-up theta
close all
font_size = 10;
fig_n = 2;
fig = figure(fig_n);
rowN = 3; colN = 2;
set(gcf,'units','centimeters','position',[0 0 colN*10 rowN*7]);

col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = 0.02; % gap width between subplots
gap_h = 0.04; % gap height between subplots
marg_b = 0.04; % bottom margin
marg_t = 0.04; % top margin
marg_l = 0.04; % left margin
marg_r = 0.02; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length
h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

% 1)
sp = 1;
axes(h_axes_sp(sp))
depth_lvl = 0;
depth_lvl_ind = depth == depth_lvl;
data = aus8_thet.Jan(:,:,depth_lvl_ind);
data_min = 7;
data_max = 21;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jan $\Theta$ $z=' num2str(depth_lvl) '$ ($^{\circ}C$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 2)
sp = 2;
axes(h_axes_sp(sp))
depth_lvl = 0;
depth_lvl_ind = depth == depth_lvl;
data = aus8_thet.Jul(:,:,depth_lvl_ind);
data_min = 7;
data_max = 21;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jul $\Theta$ $z=' num2str(depth_lvl) '$ ($^{\circ}C$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 3)
sp = 3;
axes(h_axes_sp(sp))
depth_lvl = -400;
depth_lvl_ind = depth == depth_lvl;
data = aus8_thet.Jan(:,:,depth_lvl_ind);
data_min = 7;
data_max = 11;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jan $\Theta$ $z=' num2str(depth_lvl) '$ ($^{\circ}C$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 4)
sp = 4;
axes(h_axes_sp(sp))
depth_lvl = -400;
depth_lvl_ind = depth == depth_lvl;
data = aus8_thet.Jul(:,:,depth_lvl_ind);
data_min = 7;
data_max = 11;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jul $\Theta$ $z=' num2str(depth_lvl) '$ ($^{\circ}C$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 5)
sp = 5;
axes(h_axes_sp(sp))
depth_lvl = -1000;
depth_lvl_ind = depth == depth_lvl;
data = aus8_thet.Jan(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jan $\Theta$ $z=' num2str(depth_lvl) '$ ($^{\circ}C$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 6)
sp = 6;
axes(h_axes_sp(sp))
depth_lvl = -1000;
depth_lvl_ind = depth == depth_lvl;
data = aus8_thet.Jul(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jul $\Theta$ $z=' num2str(depth_lvl) '$ ($^{\circ}C$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% Save
outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
export_fig(fig, ...
    [figures_path mfilename '/' scriptname(1:3) ...
    '_fig' num2str(fig_n) '_'], ...
    '-m3')
close


%% figure set-up asal
close all
font_size = 10;
fig_n = 3;
fig = figure(fig_n);
rowN = 3; colN = 2;
set(gcf,'units','centimeters','position',[0 0 colN*10 rowN*7]);

col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = 0.02; % gap width between subplots
gap_h = 0.04; % gap height between subplots
marg_b = 0.04; % bottom margin
marg_t = 0.04; % top margin
marg_l = 0.04; % left margin
marg_r = 0.02; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length
h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

% 1)
sp = 1;
axes(h_axes_sp(sp))
depth_lvl = 0;
depth_lvl_ind = depth == depth_lvl;
data = aus8_asal.Jan(:,:,depth_lvl_ind);
data_min = 34;
data_max = 37;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('PuOr8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jan $S_{A}$ $z=' num2str(depth_lvl) '$ ($g/kg$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 2)
sp = 2;
axes(h_axes_sp(sp))
depth_lvl = 0;
depth_lvl_ind = depth == depth_lvl;
data = aus8_asal.Jul(:,:,depth_lvl_ind);
data_min = 34;
data_max = 37;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('PuOr8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jul $S_{A}$ $z=' num2str(depth_lvl) '$ ($g/kg$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 3)
sp = 3;
axes(h_axes_sp(sp))
depth_lvl = -400;
depth_lvl_ind = depth == depth_lvl;
data = aus8_asal.Jan(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = 35.2;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('PuOr8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jan $S_{A}$ $z=' num2str(depth_lvl) '$ ($g/kg$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 4)
sp = 4;
axes(h_axes_sp(sp))
depth_lvl = -400;
depth_lvl_ind = depth == depth_lvl;
data = aus8_asal.Jul(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = 35.2;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('PuOr8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jul $S_{A}$ $z=' num2str(depth_lvl) '$ ($g/kg$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 5)
sp = 5;
axes(h_axes_sp(sp))
depth_lvl = -1000;
depth_lvl_ind = depth == depth_lvl;
data = aus8_asal.Jan(:,:,depth_lvl_ind);
data_min = 34.5;
data_max = 34.62;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('PuOr8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jan $S_{A}$ $z=' num2str(depth_lvl) '$ ($g/kg$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 6)
sp = 6;
axes(h_axes_sp(sp))
depth_lvl = -1000;
depth_lvl_ind = depth == depth_lvl;
data = aus8_asal.Jul(:,:,depth_lvl_ind);
data_min = 34.5;
data_max = 34.62;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('PuOr8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jul $S_{A}$ $z=' num2str(depth_lvl) '$ ($g/kg$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% Save
outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
export_fig(fig, ...
    [figures_path mfilename '/' scriptname(1:3) ...
    '_fig' num2str(fig_n) '_'], ...
    '-m3')
close


%% figure set-up rho
close all
font_size = 10;
fig_n = 3;
fig = figure(fig_n);
rowN = 3; colN = 2;
set(gcf,'units','centimeters','position',[0 0 colN*10 rowN*7]);

col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = 0.02; % gap width between subplots
gap_h = 0.04; % gap height between subplots
marg_b = 0.04; % bottom margin
marg_t = 0.04; % top margin
marg_l = 0.04; % left margin
marg_r = 0.02; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length
h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

% 1)
sp = 1;
axes(h_axes_sp(sp))
depth_lvl = 0;
depth_lvl_ind = depth == depth_lvl;
data = aus8_rho.Jan(:,:,depth_lvl_ind);
data_min = 1025;
data_max = 1026.8;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('BrBG8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jan $\rho$ $z=' num2str(depth_lvl) '$ ($kg/m^3$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 2)
sp = 2;
axes(h_axes_sp(sp))
depth_lvl = 0;
depth_lvl_ind = depth == depth_lvl;
data = aus8_rho.Jul(:,:,depth_lvl_ind);
data_min = 1025;
data_max = 1026.8;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('BrBG8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jul $\rho$ $z=' num2str(depth_lvl) '$ ($kg/m^3$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 3)
sp = 3;
axes(h_axes_sp(sp))
depth_lvl = -400;
depth_lvl_ind = depth == depth_lvl;
data = aus8_rho.Jan(:,:,depth_lvl_ind);
data_min = 1028.5;
data_max = 1028.7;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('BrBG8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jan $\rho$ $z=' num2str(depth_lvl) '$ ($kg/m^3$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 4)
sp = 4;
axes(h_axes_sp(sp))
depth_lvl = -400;
depth_lvl_ind = depth == depth_lvl;
data = aus8_rho.Jul(:,:,depth_lvl_ind);
data_min = 1028.5;
data_max = 1028.7;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('BrBG8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jul $\rho$ $z=' num2str(depth_lvl) '$ ($kg/m^3$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 5)
sp = 5;
axes(h_axes_sp(sp))
depth_lvl = -1000;
depth_lvl_ind = depth == depth_lvl;
data = aus8_rho.Jan(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('BrBG8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jan $\rho$ $z=' num2str(depth_lvl) '$ ($kg/m^3$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% 6)
sp = 6;
axes(h_axes_sp(sp))
depth_lvl = -1000;
depth_lvl_ind = depth == depth_lvl;
data = aus8_rho.Jul(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('BrBG8', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8 Jul $\rho$ $z=' num2str(depth_lvl) '$ ($kg/m^3$)'])
grid
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')
if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end

% Save
outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
export_fig(fig, ...
    [figures_path mfilename '/' scriptname(1:3) ...
    '_fig' num2str(fig_n) '_'], ...
    '-m3')
close



%% save the structure
save([data_path 'SACS_data/aus8_asal'], 'aus8_asal')
save([data_path 'SACS_data/aus8_thet'], 'aus8_thet')
save([data_path 'SACS_data/aus8_rho'], 'aus8_rho')

disp('aus8_asal aus8_thet aus8_rho DONE')

