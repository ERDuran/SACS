%%
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_currents'])
load([data_path 'SACS_data/KDau_fcrt'])
load([data_path 'SACS_data/aus8_figures'])
load([data_path 'SACS_data/KDau_fulu'])
load([data_path 'SACS_data/KDau_fulv'])
load([data_path 'SACS_data/KDau_salt'])
load([data_path 'SACS_data/aus8_salt'])
load([data_path 'SACS_data/aus8_u_g_prime'])
load([data_path 'SACS_data/aus8_v_g_prime'])

MTH = aus8_coor.MTH;
lat = aus8_coor.lat;
lon = aus8_coor.lon;
depth_mid = aus8_coor.depth_mid;
depth = aus8_coor.depth;
depth_thkn = aus8_coor.depth_thkn;

lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;


%% 5) plot maps of U and V SBC
N = 47;

salt_3 = KDau_salt.Jan(:,:,N:N+1);
salt_3(:,:,:,2) = KDau_salt.Feb(:,:,N:N+1);
salt_3(:,:,:,3) = KDau_salt.Mar(:,:,N:N+1);
salt_mean_KDau_temp = nanmean(salt_3,4);
salt_3 = aus8_salt.Jan(:,:,N:N+1);
salt_3(:,:,:,2) = aus8_salt.Feb(:,:,N:N+1);
salt_3(:,:,:,3) = aus8_salt.Mar(:,:,N:N+1);
salt_mean_aus8_temp = nanmean(salt_3,4);
[salt_mean_KDau, salt_mean_aus8] = deal(NaN(length(lat), length(lon)), 2);
for i = 1 : length(lat)
    for j = 1 : length(lon)
        salt_mean_KDau(i,j,1) = ...
            interp1(1:2, squeeze(salt_mean_KDau_temp(i,j,:)), 1.5);
        salt_mean_aus8(i,j,1) = ...
            interp1(1:2, squeeze(salt_mean_aus8_temp(i,j,:)), 1.5);
    end
    disp(['lat = ' num2str(i)])
end

salt_3 = KDau_salt.Apr(:,:,N:N+1);
salt_3(:,:,:,2) = KDau_salt.May(:,:,N:N+1);
salt_3(:,:,:,3) = KDau_salt.Jun(:,:,N:N+1);
salt_mean_KDau_temp = nanmean(salt_3,4);
salt_3 = aus8_salt.Apr(:,:,N:N+1);
salt_3(:,:,:,2) = aus8_salt.May(:,:,N:N+1);
salt_3(:,:,:,3) = aus8_salt.Jun(:,:,N:N+1);
salt_mean_aus8_temp = nanmean(salt_3,4);
for i = 1 : length(lat)
    for j = 1 : length(lon)
        salt_mean_KDau(i,j,2) = ...
            interp1(1:2, squeeze(salt_mean_KDau_temp(i,j,:)), 1.5);
        salt_mean_aus8(i,j,2) = ...
            interp1(1:2, squeeze(salt_mean_aus8_temp(i,j,:)), 1.5);
    end
    disp(['lat = ' num2str(i)])
end


u_g_3 = aus8_u_g_prime.Jan(:,:,N);
u_g_3(:,:,2) = aus8_u_g_prime.Feb(:,:,N);
u_g_3(:,:,3) = aus8_u_g_prime.Mar(:,:,N);
u_g_mean_aus8 = nanmean(u_g_3,3);
v_g_3 = aus8_v_g_prime.Jan(:,:,N);
v_g_3(:,:,2) = aus8_v_g_prime.Feb(:,:,N);
v_g_3(:,:,3) = aus8_v_g_prime.Mar(:,:,N);
v_g_mean_aus8 = nanmean(v_g_3,3);

u_g_3 = aus8_u_g_prime.Apr(:,:,N);
u_g_3(:,:,2) = aus8_u_g_prime.May(:,:,N);
u_g_3(:,:,3) = aus8_u_g_prime.Jun(:,:,N);
u_g_mean_aus8(:,:,2) = nanmean(u_g_3,3);
v_g_3 = aus8_v_g_prime.Apr(:,:,N);
v_g_3(:,:,2) = aus8_v_g_prime.May(:,:,N);
v_g_3(:,:,3) = aus8_v_g_prime.Jun(:,:,N);
v_g_mean_aus8(:,:,2) = nanmean(v_g_3,3);

data_pcolor = {salt_mean_KDau, salt_mean_aus8};


%%
screen_ratio = 0.75;
fig_n = 1;
rowcols = [4 1];
rowcols_size = [14 4.1]/screen_ratio; % cm
margs = [0.9 0.2 1.8 0.6]/screen_ratio; % cm
gaps = [0.4 0.6]/screen_ratio; % cm
plot_cbar_gap = 1/screen_ratio;
cbar_x = rowcols_size(1);
cbar_y = 0.2/screen_ratio;

magnif = 1;
cmap1_cont = -[200 100 60 20 5 2 0.5 0]*magnif;
cmap2_cont = -fliplr(cmap1_cont);
lvl_cmap1 = length(cmap1_cont)-1;
lvl_cmap2 = length(cmap2_cont)-1;
cmap1 = flipud(othercolor('Blues8', lvl_cmap1));
cmap2 = othercolor('Reds8', lvl_cmap2);
cmap1(end,:) = [1 1 1];
cmap2(1,:) = [1 1 1];
cmaps = [cmap1; cmap2];
cmaps_cont = [cmap1_cont cmap2_cont(2:end)];
cmaps_cont_length = length(cmaps_cont);
cmaps_linspace = linspace(0,1,cmaps_cont_length);
cmaps_y_label = cmaps_cont/magnif;

% lon_min = 115; lon_max = 147; lat_min = -47; lat_max = -32;
lon_min = 111; lon_max = 151; lat_min = -47; lat_max = -33;

x_chc = {aus8_coor.lon_u, aus8_coor.lon_v};
x_ind = [1 1];
y_chc = {aus8_coor.lat_u, aus8_coor.lat_v};

data = ...
    {interp2(lon_u,lat_u,KDau_fulu.JFM(:,:,N),lon_v,lat_v)*magnif, ...
    interp2(lon_u,lat_u,KDau_fulu.AMJ(:,:,N),lon_v,lat_v)*magnif, ...
    interp2(lon_u,lat_u,u_g_mean_aus8(:,:,1),lon_v,lat_v)*magnif, ...
    interp2(lon_u,lat_u,u_g_mean_aus8(:,:,2),lon_v,lat_v)*magnif};
v_data = {KDau_fulv.JFM(:,:,N)*magnif, KDau_fulv.AMJ(:,:,N)*magnif, ...
    v_g_mean_aus8(:,:,1)*magnif, v_g_mean_aus8(:,:,2)*magnif};

minmax = [cmaps_cont(1) cmaps_cont(end)];
cmaps_custom = cmapcust(cmaps,cmaps_cont);

axis_setup = [lon_min lon_max lat_min lat_max];
x = x_chc{x_ind(1)};
y = y_chc{x_ind(1)};

lett = 'a':'z';
font_size = 8*screen_ratio;
nan_color = [0.7 0.7 0.7];
fig_color = [1 1 1];

close all
fig = figure(fig_n);
rowN = rowcols(1); colN = rowcols(2);
[rm, cm] = meshgrid(rowN:-1:1, 1:colN);
x_sp = rowcols_size(1); % x subplot length
y_sp = rowcols_size(2); % y subplot length
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = gaps(1); % gap width between subplots
gap_h = gaps(2); % gap height between subplots
marg_b = margs(3); % bottom margin
marg_t = margs(4); % top margin
marg_l = margs(1); % left margin
marg_r = margs(2); % right margin
fig_x = marg_l+colN*x_sp+gap_w*(colN-1)+marg_r;
fig_y = marg_b+rowN*y_sp+gap_h*(rowN-1)+marg_t;

desired_length = 0.05/screen_ratio; %cm
if y_sp > x_sp, long_side = y_sp; else, long_side = x_sp; end
norm_length = desired_length/long_side;
fig_tick_length = [norm_length; 0.01];
if cbar_y > cbar_x, long_side = cbar_y; else, long_side = cbar_x; end
norm_length = desired_length/long_side;
cbar_tick_length = norm_length;

set(fig,'units','centimeters','paperunits','centimeters', ...
    'inverthardcopy','off','color',fig_color,...
    'paperposition',[0 0 fig_x fig_y]*screen_ratio,...
    'position',[0 0 fig_x fig_y]);

season = {'Summer', 'Autumn'};
for sp = 1:4
%     cbar_contour = linspace(34.40,34.60,16);
    cbar_contour = linspace(34.20,34.70,21);
    
    disp(depth_mid(N))
    subplot_x = marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1);
    subplot_y = marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1);
    ax = axes('Units','centimeters', ...
        'Position',[subplot_x,subplot_y,x_sp,y_sp]);
    
    % create reference arrow
    lat_ref = -34; lon_ref = 144;
    lat_ref_ind = find(lat_u==lat_ref+1/16); lon_ref_ind = find(lon_u==lon_ref);
    
    c = 0.1;
    data{sp}(lat_ref_ind-4:lat_ref_ind,lon_ref_ind:lon_ref_ind+2) = c;
    v_data{sp}(lat_ref_ind-15:lat_ref_ind,lon_ref_ind:lon_ref_ind+32) = 0;
    
    
    cmap = colormap(ax, flipud(othercolor('Spectral11')));
    
    if sp <= 2
        CC = 1;
    else
        CC = 2;
    end
    if sp == 1 || sp == 3
        DD = 1;
    else
        DD = 2;
    end
    
    contourf(lon,lat,data_pcolor{CC}(:,:,DD),cbar_contour,'linewidth',0.1)
    axis(axis_setup)
    shading interp
    get_caxis = caxis;
    
    text(143.5, -34.5, [num2str(c*100) ' $cm/s$'], 'color', 'k', ...
        'fontsize',font_size)
    
    hold on
    if sp == 1
        s = 3;
    elseif sp == 2
        s = 2.5;
    elseif sp == 3
        s = 2.4;
    elseif sp == 4
        s = 2.2;
    end
    
    n = 6; m = 5;
    quiver(lon_v(1:n:end),lat_v(1:m:end),...
        data{sp}(1:m:end,1:n:end),...
        v_data{sp}(1:m:end,1:n:end),...
        s,'k')
    
    if sp == 1
    h_tit = title(['(' lett(sp) ') MOM01-75z summer ' ...
        ' salinity and velocity field at z = ' ...
        num2str(depth_mid(N)) ' $m$'], ...
        'horizontalalignment','left', 'fontsize',font_size);
    elseif sp == 3
    h_tit = title(['(' lett(sp) ') CARS-aus8 summer ' ...
        ' salinity and velocity field at z = ' ...
        num2str(depth_mid(N)) ' $m$'], ...
        'horizontalalignment','left', 'fontsize',font_size);
    elseif sp == 2
    h_tit = title(['(' lett(sp) ') MOM01-75z autumn ' ...
        ' salinity and velocity field at z = ' ...
        num2str(depth_mid(N)) ' $m$'], ...
        'horizontalalignment','left', 'fontsize',font_size);   
    elseif sp == 4
    h_tit = title(['(' lett(sp) ') CARS-aus8 autumn ' ...
        ' salinity and velocity field at z = ' ...
        num2str(depth_mid(N)) ' $m$'], ...
        'horizontalalignment','left', 'fontsize',font_size);   
    end
    
    h_tit.Position(1) = axis_setup(1);
    grid
    set(ax,'layer','top','color',nan_color,...
        'fontsize',font_size,'tickdir','out', ...
        'ticklength',fig_tick_length, ...
        'xtick', lon_min:2:lon_max, 'ytick', -50:2:-30)
    if row_ind(sp) ~= rowN
        set(gca,'xticklabel','')
    else
        xlabel('Longitude')
    end
    if col_ind(sp) ~= 1
        set(gca,'yticklabel','')
    else
        ylabel('Latitude')
    end
    
    if sp == 4
        ax = axes('visible', 'off');
        
        cmap = colormap(ax, ...
            flipud(othercolor('Spectral11', length(cbar_contour)-1)));
        colormap(ax, cmap);
        cbar = colorbar('horizontal');
%         set(cbar, ...
%             'YAxisLocation','right',...
%             'ytick',linspace(0,1,length(cbar_contour)), ...
%             'yticklabel',round(cbar_contour*100)/100,...
%             'fontsize',font_size,'ticklength',cbar_tick_length)
        set(cbar, ...
            'YAxisLocation','right',...
            'ytick',linspace(0,1,length(cbar_contour)), ...
            'yticklabel',cbar_contour,...
            'fontsize',font_size,'ticklength',cbar_tick_length)
        set(cbar,'units','centimeters','position', [...
            (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)), ...
            (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)-plot_cbar_gap), ...
            cbar_x, ...
            cbar_y]);
        set(get(cbar,'xlabel'),'String','Salinity ($psu$)', ...
            'fontsize',font_size)
        cbar.Label.Interpreter = 'latex';
        
        pointycbar(cbar)
    end
end

outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
print(fig, ...
    [figures_path mfilename '/' scriptname(1:3) ...
    '_fig' num2str(fig_n) '_'], ...
    '-dpng', '-r300')
print(fig, ...
    ['~/Duran2017/SACS/10319442jbhpxfsdfvwy/' scriptname(1:3) ...
    '_fig' num2str(fig_n) '_'], ...
    '-dpng', '-r300')
close


%%
play(bird_f_path,[1 (get(bird_f_path, 'SampleRate')*3)]);

