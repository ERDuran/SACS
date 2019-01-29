%%
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_currents'])
load([data_path 'SACS_data/aus8_U_prime'])
load([data_path 'SACS_data/aus8_V_prime'])

a = aus8_coor.a;
pi180 = aus8_coor.pi180;
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;
lat_v_SBC_north = aus8_currents.lat_v_SBC_north;
lat_v_SBC_south = aus8_currents.lat_v_SBC_south;
lat_v_DRC_north = aus8_currents.lat_v_DRC_north;
lat_v_DRC_south = aus8_currents.lat_v_DRC_south;
lon_u_ALLC = aus8_currents.lon_u_ALLC;
Months = aus8_coor.Months;
Months{13} = 'mean';


%% Prep
DRC_Ut_up_temp = NaN(12,361);

for t = 1 : 12
    %
    U_prime_up = aus8_currents.ztop_to_zmid.U_prime.(Months{t});
    V_prime_up = aus8_currents.ztop_to_zmid.V_prime.(Months{t});
    U_g_prime_dw = aus8_currents.zmid_to_zbot.U_g_prime.(Months{t});
    V_g_prime_dw = aus8_currents.zmid_to_zbot.V_g_prime.(Months{t});
    
    %
    dx_v = NaN(length(lat_v), length(lon_u)-1);
    for ii = 1 : length(lat_v)
        dx_v(ii,:) = a * cos(lat_v(ii) * pi180) * ...
            (lon_u(2:end) - lon_u(1:end-1)) * pi180;
    end
    dy_raw = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
    dy_u = repmat(dy_raw, [1 length(lon_u)]);
    
    %
    U_prime_up(isnan(U_prime_up)) = 0;
    V_prime_up(isnan(V_prime_up)) = 0;
    U_g_prime_dw(isnan(U_g_prime_dw)) = 0;
    V_g_prime_dw(isnan(V_g_prime_dw)) = 0;
    Ut_prime_up = U_prime_up .* dy_u;
    Vt_prime_up = V_prime_up .* dx_v;
    Ut_g_prime_dw = U_g_prime_dw .* dy_u;
    Vt_g_prime_dw = V_g_prime_dw .* dx_v;
    
    %
    lon_u_ALLC_repelem = [...
        lon_u_ALLC(:,1), ...
        repelem(lon_u_ALLC(:,2:end-1), 1, 2), ...
        lon_u_ALLC(:,end)];
    lat_v_SBC_north_repelem = ...
        repelem(lat_v_SBC_north, 1, 2);
    lat_v_SBC_south_repelem = ...
        repelem(lat_v_SBC_south, 1, 2);
    lat_v_DRC_north_repelem = ...
        repelem(lat_v_DRC_north, 1, 2);
    lat_v_DRC_south_repelem = ...
        repelem(lat_v_DRC_south, 1, 2);
    
    %
    dx_u = NaN(length(lat_u), length(lon_u)-1);
    for ii = 1 : length(lat_u)
        dx_u(ii,:) = a * cos(lat_u(ii) * pi180) * ...
            (lon_u(2:end) - lon_u(1:end-1)) * pi180;
    end
    dy_raw = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
    dy_v = repmat(dy_raw, [1 length(lon_v)]);
    
    %
    aus8_U_prime.(Months{t})(isnan(aus8_U_prime.(Months{t}))) = 0;
    aus8_V_prime.(Months{t})(isnan(aus8_V_prime.(Months{t}))) = 0;
    du = aus8_U_prime.(Months{t})(:,2:end) - ...
        aus8_U_prime.(Months{t})(:,1:end-1);
    dudx = du ./ dx_u;
    lat_v_repmat = repmat(lat_v,1,length(lon_v));
    dv = aus8_V_prime.(Months{t})(1:end-1,:).* ...
        cos(lat_v_repmat(1:end-1,:) * pi180) - ...
        aus8_V_prime.(Months{t})(2:end,:).* ...
        cos(lat_v_repmat(2:end,:) * pi180);
    lat_u_repmat = repmat(lat_u,1,length(lon_v));
    dvdy = 1./cos(lat_u_repmat * pi180).*dv./dy_v;
    div_UV_prime.(Months{t}) = dudx + dvdy;
    
    % DRC Ut up
    Ut_prime_up_DRC = zeros(size(Ut_prime_up));
    
    % First
    lat_north_ind_1 = ...
        find(ismember(lat_v, lat_v_DRC_north_repelem(1)));
    lat_south_ind_1 = ...
        find(ismember(lat_v, lat_v_DRC_south_repelem(1))) -1;
    lat_vec_now = lat_north_ind_1:lat_south_ind_1;
    lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(1));
    Ut_prime_up_DRC(lat_vec_now, ...
        lon_u_ALLC_repelem_ind) = ...
        Ut_prime_up(lat_vec_now,...
        lon_u_ALLC_repelem_ind);
    % Last
    lat_north_ind_1 = ...
        find(ismember(lat_v, lat_v_DRC_north_repelem(end)));
    lat_south_ind_1 = ...
        find(ismember(lat_v, lat_v_DRC_south_repelem(end))) -1;
    lat_vec_now = lat_north_ind_1:lat_south_ind_1;
    lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(end));
    Ut_prime_up_DRC(lat_vec_now, ...
        lon_u_ALLC_repelem_ind) = ...
        Ut_prime_up(lat_vec_now,...
        lon_u_ALLC_repelem_ind);
    
    jj_count = 0;
    for jj = length(lon_u_ALLC_repelem(1:end-1)) : -2 : 3
        jj_count = jj_count + 1;
        
        lat_north_ind_1 = ...
            find(ismember(lat_v, lat_v_DRC_north_repelem(jj)));
        lat_north_ind_2 = ...
            find(ismember(lat_v, lat_v_DRC_north_repelem(jj-1)));
        lat_north_12 = [lat_north_ind_1, lat_north_ind_2];
        lat_south_ind_1 = ...
            find(ismember(lat_v, lat_v_DRC_south_repelem(jj))) -1;
        lat_south_ind_2 = ...
            find(ismember(lat_v, lat_v_DRC_south_repelem(jj-1))) -1;
        lat_south_12 = [lat_south_ind_1, lat_south_ind_2];
        
        if lat_north_12(1) == lat_north_12(2)
            north_12_ind = 1;
        else
            north_12_ind = find(lat_north_12 == max(lat_north_12));
        end
        if lat_south_12(1) == lat_south_12(2)
            south_12_ind = 1;
        else
            south_12_ind = find(lat_south_12 == min(lat_south_12));
        end
        
        lat_vec_now = ...
            lat_north_12(north_12_ind):lat_south_12(south_12_ind);
        lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(jj));
        Ut_prime_up_DRC(lat_vec_now, ...
            lon_u_ALLC_repelem_ind) = ...
            Ut_prime_up(lat_vec_now,...
            lon_u_ALLC_repelem_ind);
    end
    
    DRC_Ut_mcps_up = sum(Ut_prime_up_DRC, 1);
    DRC_Ut_up_temp(t,:) = DRC_Ut_mcps_up * 10^-6;
end


%% 5) plot maps of U and V SBC
screen_ratio = 0.75;
fig_n = 1;
rowcols = [1 1];
rowcols_size = [14 6]/screen_ratio/2; % cm
margs = [1.0 0.5 2.9 0.7]/screen_ratio; % cm
gaps = [0.4 0.7]/screen_ratio; % cm

magnif = 10;
DRC_Ut_up = DRC_Ut_up_temp*magnif;
cmap1_cont = -[4 3 2 1 0.5 0.1 0]*magnif;
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
lon_min = 110; lon_max = 152; lat_min = -48; lat_max = -31;

x_chc = {aus8_coor.lon_u, aus8_coor.lon_v, aus8_coor.lon_u, aus8_coor.lon_v};
x_ind = [1 1 1 1];
y_chc = {aus8_coor.lat_u, aus8_coor.lat_v, aus8_coor.lat_u, aus8_coor.lat_v};

for sp = 1 : rowcols(1)*rowcols(2)
    minmax{sp} = [cmaps_cont(1) cmaps_cont(end)];
    cmaps_custom{sp} = cmapcust(cmaps,cmaps_cont);
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

plot_cbar_gap = -2.1/screen_ratio;
cbar_x = x_sp;
cbar_y = 0.2/screen_ratio;

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
    
    colormap(ax, cmaps_custom{sp});
    pcolor(lon_u, 1:12, DRC_Ut_up)
    axis([115 147 1 12])
    shading interp
    caxis([minmax{sp}(1) minmax{sp}(2)]);
    
    h_tit = title(['CARS-aus8 monthly $\mathcal{U}_{FC}$ in the upper layer from the surface to 250 $m$'], ...
        'horizontalalignment','left', 'fontsize',font_size, ...
        'Interpreter','latex');
    h_tit.Position(1) = 115;
    grid
    set(ax,'layer','top','color',nan_color,...
        'fontsize',font_size,'tickdir','out', ...
        'ticklength',fig_tick_length, ...
        'xtick', 115:2:147, 'ytick', 1:12, 'yticklabel', Months)
    if row_ind(sp) ~= rowN
        set(gca,'xticklabel','')
    else
        xlabel('Longitude', 'fontsize', font_size)
    end
    
    ylabel('Month', 'fontsize', font_size)
    
    small_arr = 3;
    % Cape Leeuwin
    hh = arrow([115,-0.25], [115,-1.75], small_arr, ...
        'Facecolor', 'k', 'edgecolor', 'k');
    text(115,-2.2,'Cape Leeuwin', 'fontsize', font_size, ...
        'color', 'k','Rotation',-45)
    set(get(get(hh,'Annotation'),'LegendInformation'), ...
        'IconDisplayStyle','off');
    % Albany
    hh = arrow([118,-0.25], [118,-1.75], small_arr, ...
        'Facecolor', 'k', 'edgecolor', 'k');
    text(118,-2.2,'Albany', 'fontsize', font_size, ...
        'color', 'k','Rotation',-45)
    set(get(get(hh,'Annotation'),'LegendInformation'), ...
        'IconDisplayStyle','off');
    % Cape Pasley
    hh = arrow([123.4,-0.25], [123.4,-1.75], small_arr, ...
        'Facecolor', 'k', 'edgecolor', 'k');
    text(123.4,-2.2,'Cape Pasley', 'fontsize', font_size, ...
        'color', 'k','Rotation',-45)
    set(get(get(hh,'Annotation'),'LegendInformation'), ...
        'IconDisplayStyle','off');
    % West GAB
    hh = arrow([128,-0.25], [128,-1.75], small_arr, ...
        'Facecolor', 'k', 'edgecolor', 'k');
    text(128,-2.2,'West GAB', 'fontsize', font_size, ...
        'color', 'k','Rotation',-45)
    set(get(get(hh,'Annotation'),'LegendInformation'), ...
        'IconDisplayStyle','off');
    % East GAB
    hh = arrow([132,-0.25], [132,-1.75], small_arr, ...
        'Facecolor', 'k', 'edgecolor', 'k');
    text(132,-2.2,'East GAB', 'fontsize', font_size, ...
        'color', 'k','Rotation',-45)
    set(get(get(hh,'Annotation'),'LegendInformation'), ...
        'IconDisplayStyle','off');
    % Cape Carnot
    hh = arrow([135.7,-0.25], [135.7,-1.75], small_arr, ...
        'Facecolor', 'k', 'edgecolor', 'k');
    text(135.7,-2.2,'Cape Carnot', 'fontsize', font_size, ...
        'color', 'k','Rotation',-45)
    set(get(get(hh,'Annotation'),'LegendInformation'), ...
        'IconDisplayStyle','off');
    % Portland
    hh = arrow([141.3,-0.25], [141.3,-1.75], small_arr, ...
        'Facecolor', 'k', 'edgecolor', 'k');
    text(141.3,-2.2,'Portland', 'fontsize', font_size, ...
        'color', 'k','Rotation',-45)
    set(get(get(hh,'Annotation'),'LegendInformation'), ...
        'IconDisplayStyle','off');
    %  King Is.
    hh = arrow([144,-0.25], [144,-1.75], small_arr, ...
        'Facecolor', 'k', 'edgecolor', 'k');
    text(144,-2.2,' King Island', 'fontsize', font_size, ...
        'color', 'k','Rotation',-45)
    set(get(get(hh,'Annotation'),'LegendInformation'), ...
        'IconDisplayStyle','off');
    % South East Cape
    hh = arrow([146.7,-0.25], [146.7,-1.75], small_arr, ...
        'Facecolor', 'k', 'edgecolor', 'k');
    text(146.7,-2.2,'SE C.', 'fontsize', font_size, ...
        'color', 'k','Rotation',-45)
    set(get(get(hh,'Annotation'),'LegendInformation'), ...
        'IconDisplayStyle','off');
    
    if sp == 1
        ax = axes('visible', 'off');
        colormap(ax, cmaps);
        cbar = colorbar('horizontal');
        set(cbar,'xtick',cmaps_linspace, ...
            'XTickLabel',cmaps_y_label,...
            'fontsize',font_size,'ticklength',cbar_tick_length)
        set(cbar,'units','centimeters','position', [...
            marg_l, ...
            marg_b+plot_cbar_gap, ...
            cbar_x, ...
            cbar_y]);
        
        set(get(cbar,'xlabel'), ...
            'String','Transport ($Sv$)', ...
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

