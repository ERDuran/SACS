%% make monthly variations
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8'])
load([data_path 'SACS_data/aus8_coor'])
lat = aus8_coor.lat;
lon = aus8_coor.lon;
depth = aus8_coor.depth;


%% 
time_step = [...
    15, 46, 74, 105, 135, 166, ...
    196, 227, 258, 288, 319, 349];

month_names = {...
    'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

% time cycle setup
time_serie = 2 * pi / 365 * time_step;

for ll = 1 : 12
    temp_monthly = NaN(size(aus8.temp.mean));
    salt_monthly = NaN(size(aus8.salt.mean));
    
    for kk = 1 : length(depth)
        if depth(kk) >= aus8.depth_semiann(end)
            temp_monthly(:,:,kk) = ...
                ... % mean component
                aus8.temp.mean(:,:,kk) + ...
                ... % annual component
                aus8.temp.an_cos(:,:,kk)*cos(time_serie(ll)) + ...
                aus8.temp.an_sin(:,:,kk)*sin(time_serie(ll)) + ...
                ... % semi-annual component
                aus8.temp.sa_cos(:,:,kk)*cos(2*time_serie(ll)) + ...
                aus8.temp.sa_sin(:,:,kk)*sin(2*time_serie(ll));
            
            salt_monthly(:,:,kk) = ...
                aus8.salt.mean(:,:,kk) + ...
                aus8.salt.an_cos(:,:,kk)*cos(time_serie(ll)) + ...
                aus8.salt.an_sin(:,:,kk)*sin(time_serie(ll)) + ...
                aus8.salt.sa_cos(:,:,kk)*cos(2*time_serie(ll)) + ...
                aus8.salt.sa_sin(:,:,kk)*sin(2*time_serie(ll));
            
            % depths above the annual component max depth
        elseif depth(kk) >= aus8.depth_ann(end)
            temp_monthly(:,:,kk) = ...
                ... % mean component
                aus8.temp.mean(:,:,kk) + ...
                ... % annual component
                aus8.temp.an_cos(:,:,kk)*cos(time_serie(ll)) + ...
                aus8.temp.an_sin(:,:,kk)*sin(time_serie(ll));
            
            salt_monthly(:,:,kk) = ...
                aus8.salt.mean(:,:,kk) + ...
                aus8.salt.an_cos(:,:,kk)*cos(time_serie(ll)) + ...
                aus8.salt.an_sin(:,:,kk)*sin(time_serie(ll));
            
            % anything below the max depth of the annual component
            % is the simply the mean
        else
            temp_monthly(:,:,kk) = ...
                ... % mean component
                aus8.temp.mean(:,:,kk);
            
            salt_monthly(:,:,kk) = ...
                aus8.salt.mean(:,:,kk);
        end
    end
    
    aus8_monthly.temp.(month_names{ll}) = temp_monthly;
    aus8_monthly.salt.(month_names{ll}) = salt_monthly;
    disp([month_names{ll} ' OK!'])
end


%% figure set-up
close all
font_size = 8;
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
data = aus8_monthly.temp.Jan(:,:,depth_lvl_ind);
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
data = aus8_monthly.temp.Jul(:,:,depth_lvl_ind);
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
data = aus8_monthly.temp.Jan(:,:,depth_lvl_ind);
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
data = aus8_monthly.temp.Jul(:,:,depth_lvl_ind);
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
data = aus8_monthly.temp.Jan(:,:,depth_lvl_ind);
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
data = aus8_monthly.temp.Jul(:,:,depth_lvl_ind);
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
fig_n = 2;
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
data = aus8_monthly.salt.Jan(:,:,depth_lvl_ind);
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
data = aus8_monthly.salt.Jul(:,:,depth_lvl_ind);
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
data = aus8_monthly.salt.Jan(:,:,depth_lvl_ind);
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
data = aus8_monthly.salt.Jul(:,:,depth_lvl_ind);
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
data = aus8_monthly.salt.Jan(:,:,depth_lvl_ind);
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
data = aus8_monthly.salt.Jul(:,:,depth_lvl_ind);
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


%% save aus8_coast
save([data_path 'SACS_data/aus8_monthly'], 'aus8_monthly')
disp('aus8_monthly DONE')

