%% calc t s and rho
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8'])
load([data_path 'SACS_data/aus8_monthly'])
load([data_path 'SACS_data/ETOPO5'])

lat = aus8_coor.lat;
lon = aus8_coor.lon;
depth = aus8_coor.depth;
topog = ETOPO5.topo_sm_interp;


%% increase vertical resolution
aus8_coor.bottom_depth = NaN([length(lat), length(lon)]);

aus8_temp.mean = aus8.temp.mean;
aus8_salt.mean = aus8.salt.mean;

Months = {...
    'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
aus8_coor.Months = Months;
for t = 1 : 12
    aus8_temp.(Months{t}) = aus8_monthly.temp.(Months{t});
    aus8_salt.(Months{t}) = aus8_monthly.salt.(Months{t});
end

for m = 1 : length(lat)
    for n = 1 : length(lon)
        last_finite = ...
            find(isfinite(squeeze(aus8.temp.mean(m,n,:))), ...
            1, 'last');
        % apply ETOPO5 mask
        if ~isempty(last_finite)
            aus8_coor.bottom_depth(m,n) = depth(last_finite);
            if topog(m,n) > aus8_coor.bottom_depth(m,n)
                land_idx = depth <= topog(m,n);
                aus8_temp.mean(m,n,land_idx) = NaN;
                aus8_salt.mean(m,n,land_idx) = NaN;
            end
        end
        
        above_first_NaN = ...
            find(isnan(squeeze(aus8.temp.mean(m,n,:))), ...
            1, 'first') -1;
        % interpolate cavities
        if above_first_NaN ~= last_finite
            aus8_temp.mean(m,n,:) = ...
                fixgaps(squeeze(aus8.temp.mean(m,n,:)));
            aus8_salt.mean(m,n,:) = ...
                fixgaps(squeeze(aus8.salt.mean(m,n,:)));
        end
        
        for t = 1 : 12
            if ~isempty(last_finite)
                if topog(m,n) > aus8_coor.bottom_depth(m,n)
                    aus8_temp.(Months{t})(m,n,land_idx) = NaN;
                    aus8_salt.(Months{t})(m,n,land_idx) = NaN;
                end
            end
            if above_first_NaN ~= last_finite
                aus8_temp.(Months{t})(m,n,:) = ...
                    fixgaps(squeeze(...
                    aus8_monthly.temp.(Months{t})(m,n,:)));
                aus8_salt.(Months{t})(m,n,:) = ...
                    fixgaps(squeeze(...
                    aus8_monthly.salt.(Months{t})(m,n,:)));
            end
        end
    end
    fprintf('lat = %5.3f\n', lat(m))
end


%% figure set-up
close all
font_size = 10;
fig_n = 1;
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
data = aus8_temp.mean(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu11', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.temp.mean z=' num2str(depth_lvl) ' (degc)'])
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
data = aus8_salt.mean(:,:,depth_lvl_ind);
data_min = 33.5;
data_max = 37;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('YlGnBu9', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.salt.mean z=' num2str(depth_lvl) ' (psu)'])
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
data = aus8_temp.mean(:,:,depth_lvl_ind);
data_min = 6;
data_max = 12;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu11', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.temp.mean z=' num2str(depth_lvl) ' (degc)'])
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
data = aus8_salt.mean(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = 35.2;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('YlGnBu9', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.salt.mean z=' num2str(depth_lvl) ' (psu)'])
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
data = aus8_temp.mean(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu11', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.temp.mean z=' num2str(depth_lvl) ' (degC)'])
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
data = aus8_salt.mean(:,:,depth_lvl_ind);
data_min = 34.25;
data_max = 34.55;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('YlGnBu9', levels)));
pcolor(...
    aus8_coor.lon, ...
    aus8_coor.lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.salt.mean z=' num2str(depth_lvl) ' (psu)'])
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


%% figure set-up
close all
font_size = 8;
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
data = aus8_temp.Jan(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu11', levels)));
pcolor(...
    lon, ...
    lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.temp.Jan z=' num2str(depth_lvl) ' (degC)'])
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
data = aus8_temp.Jul(:,:,depth_lvl_ind);
% data_min = min(min(data));
% data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu11', levels)));
pcolor(...
    lon, ...
    lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.temp.Jul z=' num2str(depth_lvl) ' (degC)'])
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
data = aus8_temp.Jan(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu11', levels)));
pcolor(...
    lon, ...
    lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.temp.Jan z=' num2str(depth_lvl) ' (degC)'])
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
data = aus8_temp.Jul(:,:,depth_lvl_ind);
% data_min = min(min(data));
% data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu11', levels)));
pcolor(...
    lon, ...
    lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.temp.Jul z=' num2str(depth_lvl) ' (degC)'])
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
data = aus8_temp.Jan(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu11', levels)));
pcolor(...
    lon, ...
    lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.temp.Jan z=' num2str(depth_lvl) ' (degC)'])
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
data = aus8_temp.Jul(:,:,depth_lvl_ind);
% data_min = min(min(data));
% data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('RdYlBu11', levels)));
pcolor(...
    lon, ...
    lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.temp.Jul z=' num2str(depth_lvl) ' (degC)'])
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


%% figure set-up
close all
font_size = 10;
fig_n = 3;
fig = figure(fig_n);
rowN = 3; colN = 2;
set(gcf,'units','centimeters','position',[0 0 colN*10 rowN*7]);

col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = 0.04; % gap width between subplots
gap_h = 0.04; % gap height between subplots
marg_b = 0.04; % bottom margin
marg_t = 0.04; % top margin
marg_l = 0.04; % left margin
marg_r = 0.04; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length
h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

% 1)
sp = 1;
axes(h_axes_sp(sp))
depth_lvl = 0;
depth_lvl_ind = depth == depth_lvl;
data = aus8_salt.Jan(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = 37;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('YlGnBu9', levels)));
pcolor(...
    lon, ...
    lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.salt.Jan z=' num2str(depth_lvl) ' (psu)'])
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
data = aus8_salt.Jul(:,:,depth_lvl_ind);
% data_min = min(min(data));
% data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('YlGnBu9', levels)));
pcolor(...
    lon, ...
    lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.salt.Jul z=' num2str(depth_lvl) ' (psu)'])
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
data = aus8_salt.Jan(:,:,depth_lvl_ind);
data_min = 34.4;
data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('YlGnBu9', levels)));
pcolor(...
    lon, ...
    lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.salt.Jan z=' num2str(depth_lvl) ' (psu)'])
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
data = aus8_salt.Jul(:,:,depth_lvl_ind);
% data_min = min(min(data));
% data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('YlGnBu9', levels)));
pcolor(...
    lon, ...
    lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.salt.Jul z=' num2str(depth_lvl) ' (psu)'])
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
data = aus8_salt.Jan(:,:,depth_lvl_ind);
data_min = min(min(data));
data_max = 34.5;
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('YlGnBu9', levels)));
pcolor(...
    lon, ...
    lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.salt.Jan z=' num2str(depth_lvl) ' (psu)'])
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
data = aus8_salt.Jul(:,:,depth_lvl_ind);
% data_min = min(min(data));
% data_max = max(max(data));
levels = 12;
colormap(h_axes_sp(sp), flipud(othercolor('YlGnBu9', levels)));
pcolor(...
    lon, ...
    lat, ...
    data)
shading interp
caxis([data_min data_max]);
cbar = colorbar;
title(['aus8.salt.Jul z=' num2str(depth_lvl) ' (psu)'])
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


%% save updated a8c
save([data_path 'SACS_data/aus8_coor'], 'aus8_coor')
save([data_path 'SACS_data/aus8_temp'], 'aus8_temp')
save([data_path 'SACS_data/aus8_salt'], 'aus8_salt')
disp('aus8_coor aus8_temp aus8_salt DONE')

