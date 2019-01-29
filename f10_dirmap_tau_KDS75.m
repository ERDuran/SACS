%%
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/aus8_coor'], 'aus8_coor')
load([data_path 'SACS_data/KDau_tau'], 'KDau_tau')
load([data_path 'SACS_data/KDau_fcrt'])
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
Seasons = {'Summer', 'Autumn', 'Winter', 'Spring'};


%%
tau_x.JFM = ...
    (KDau_tau.tau_x.Jan + KDau_tau.tau_x.Feb + KDau_tau.tau_x.Mar)/3;
tau_x.AMJ = ...
    (KDau_tau.tau_x.Apr + KDau_tau.tau_x.May + KDau_tau.tau_x.Jun)/3;
tau_x.JAS = ...
    (KDau_tau.tau_x.Jul + KDau_tau.tau_x.Aug + KDau_tau.tau_x.Sep)/3;
tau_x.OND = ...
    (KDau_tau.tau_x.Oct + KDau_tau.tau_x.Nov + KDau_tau.tau_x.Dec)/3;

%
tau_y.JFM = ...
    (KDau_tau.tau_y.Jan + KDau_tau.tau_y.Feb + KDau_tau.tau_y.Mar)/3;
tau_y.AMJ = ...
    (KDau_tau.tau_y.Apr + KDau_tau.tau_y.May + KDau_tau.tau_y.Jun)/3;
tau_y.JAS = ...
    (KDau_tau.tau_y.Jul + KDau_tau.tau_y.Aug + KDau_tau.tau_y.Sep)/3;
tau_y.OND = ...
    (KDau_tau.tau_y.Oct + KDau_tau.tau_y.Nov + KDau_tau.tau_y.Dec)/3;

%
for t = 1 : 4
    tau_x_interp2.(MTH{t}) = interp2(lon_v, lat_v, tau_x.(MTH{t}), ...
        lon_v,lat_u);
    tau_y_interp2.(MTH{t}) = interp2(lon_u, lat_u, tau_y.(MTH{t}), ...
        lon_v,lat_u);
    speed_interp2.(MTH{t}) = ...
        sqrt(tau_x_interp2.(MTH{t}).^2 + tau_y_interp2.(MTH{t}).^2);
end


%% 5) plot maps of U and V SBC
screen_ratio = 0.75;
fig_n = 1;
rowcols = [4 1];
rowcols_size = [14 6]/screen_ratio/2; % cm
margs = [0.9 0.3 1.8 0.6]/screen_ratio; % cm
gaps = [0.4 0.8]/screen_ratio; % cm
plot_cbar_gap = 1/screen_ratio;
cbar_x = rowcols_size(1);
cbar_y = 0.2/screen_ratio;

magnif = 1000;
cmap1_cont = [0 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225]*magnif;
lvl_cmap1 = length(cmap1_cont)-1;
cmap1 = othercolor('Purples8', lvl_cmap1);
% cmap1(end,:) = [1 1 1];
cmaps_cont_length = length(cmap1_cont);
cmaps_linspace = linspace(0,1,cmaps_cont_length);
cmaps_y_label = cmap1_cont/magnif;

% lon_min = 115; lon_max = 147; lat_min = -47; lat_max = -32;
lon_min = 110; lon_max = 152; lat_min = -48; lat_max = -31;

for t = 1 : 4
    u_data{t} = tau_x_interp2.(MTH{t})*magnif;
    v_data{t} = tau_y_interp2.(MTH{t})*magnif;
    data{t} = speed_interp2.(MTH{t})*magnif;
end

for sp = 1 : rowcols(1)*rowcols(2)
    minmax{sp} = [cmap1_cont(1) cmap1_cont(end)];
    cmaps_custom{sp} = cmapcust(cmap1,cmap1_cont);
    
    axis_setup{sp} = [lon_min lon_max lat_min lat_max];

    x{sp} = lon_v;
    y{sp} = lat_u;
end

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

for sp = 1 : rowN*colN
    subplot_x = marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1);
    subplot_y = marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1);
    ax = axes('Units','centimeters', ...
        'Position',[subplot_x,subplot_y,x_sp,y_sp]);
    
    % create reference arrow
    lat_ref = -34.5; lon_ref = 142.5;
    lat_ref_ind = find(lat_u==lat_ref+1/16); 
    lon_ref_ind = find(lon_u==lon_ref);
    
    data{sp}(lat_ref_ind-30:lat_ref_ind,lon_ref_ind:lon_ref_ind+58) = 0;
    
    if sp == 1
        s = 3;
        ref_magn = 0.25;
    end
    
    u_data{sp}(25,289) = ref_magn*magnif;
    v_data{sp}(25,289) = 0;
        
    colormap(ax, cmaps_custom{sp});
    pcolor(x{sp}, y{sp}, data{sp})
    axis(axis_setup{sp})
    shading interp
    caxis([minmax{sp}(1) minmax{sp}(2)]);
    hold on
    
    s = 1.75; n = 12; m = 12;
    quiver(x{sp}(1:n:end),y{sp}(1:m:end),...
        u_data{sp}(1:m:end,1:n:end),...
        v_data{sp}(1:m:end,1:n:end),...
        s,'k')    
    
    text(lon_u(281), lat_u(18), ...
        [num2str(ref_magn) ' $N/m^{2}$'], ...
        'fontsize',font_size)
    
%     bass_strait = round(KDau_fcrt.SBC_Vtnc.(MTH{sp})(312)*100)/100;
%     
%     if bass_strait > 0
%         text(lon_u(280), lat_u(65), ...
%             [num2str(bass_strait) '$Sv$'], 'color','red',...
%             'fontsize',font_size)
%     else
%         text(lon_u(280), lat_u(65), ...
%             [num2str(bass_strait) '$Sv$'], 'color','blue',...
%             'fontsize',font_size)
%     end
    
    h_tit = title(['(' lett(sp) ') ' Seasons{sp} ' wind stress'], ...
    'horizontalalignment','left', 'fontsize',font_size);
    h_tit.Position(1) = axis_setup{sp}(1);
    grid
    set(ax,'layer','top','color',nan_color,...
        'fontsize',font_size,'tickdir','out', ...
        'ticklength',fig_tick_length, ...
        'xtick', lon_min:4:lon_max, 'ytick', lat_min:4:lat_max)
    if row_ind(sp) ~= rowN
        set(gca,'xticklabel','')
    else
        xlabel('Longitude')
    end
    ylabel('Latitude')
    if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end
    
    if sp == rowN*colN
        ax = axes('visible', 'off');
        colormap(ax, cmap1);
        cbar = colorbar('horizontal');
        set(cbar,'ytick',cmaps_linspace, ...
            'XTickLabel',cmaps_y_label,...
            'fontsize',font_size,'ticklength',cbar_tick_length)
        set(cbar,'units','centimeters','position', [...
            (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)), ...
            (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)-plot_cbar_gap), ...
            cbar_x, ...
            cbar_y]);
        set(get(cbar,'xlabel'),'String','$\tau$ ($N/m^{2}$)', ...
            'fontsize',font_size)
        cbar.Label.Interpreter = 'latex';
        
%         pointycbar(cbar)
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



