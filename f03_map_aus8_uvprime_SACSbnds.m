%% fig 1: map of SBC currents
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_currents'])


%%
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;

U_prime_ptop_to_pmid = aus8_currents.ptop_to_pmid.U_prime.mean;
V_prime_ptop_to_pmid = aus8_currents.ptop_to_pmid.V_prime.mean;
U_g_prime_pmid_to_pbot = aus8_currents.pmid_to_pbot.U_g_prime.mean;
V_g_prime_pmid_to_pbot = aus8_currents.pmid_to_pbot.V_g_prime.mean;

lon_u_ALLC_repelem = [...
    aus8_currents.lon_u_ALLC(:,1), ...
    repelem(aus8_currents.lon_u_ALLC(:,2:end-1), 1, 2), ...
    aus8_currents.lon_u_ALLC(:,end)];
lat_v_SBC_centr_repelem = ...
    repelem(aus8_currents.lat_v_SBC_centr, 1, 2);
lat_v_SBC_north_repelem = ...
    repelem(aus8_currents.lat_v_SBC_north, 1, 2);
lat_v_SBC_south_repelem = ...
    repelem(aus8_currents.lat_v_SBC_south, 1, 2);
lat_v_DRC_centr_repelem = ...
    repelem(aus8_currents.lat_v_DRC_centr, 1, 2);
lat_v_DRC_north_repelem = ...
    repelem(aus8_currents.lat_v_DRC_north, 1, 2);
lat_v_DRC_south_repelem = ...
    repelem(aus8_currents.lat_v_DRC_south, 1, 2);

z_top = aus8_currents.z_top;
z_mid = aus8_currents.z_mid;
z_bot = aus8_currents.z_bot;


%% 5) plot maps of U and V SBC
fig_n = 1;
rowcols = [2 2];
rowcols_size = [14 8]; % cm
margs = [1 2.5 1.2 1.2]; % cm
gaps = [0.75 1.25]; % cm
plot_cbar_gap = 0.5;
cbar_x = 0.3;
cbar_y = rowcols_size(2);

magnif = 10;
cmap1_cont = -[100 50 20 10 5 2 0.5 0]*magnif;
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

lon_min = 115; lon_max = 147; lat_min = -47; lat_max = -32;

x_chc = {aus8_coor.lon_u, aus8_coor.lon_v};
x_ind = [1 2 1 2];
y_chc = {aus8_coor.lat_u, aus8_coor.lat_v};

data = {U_prime_ptop_to_pmid*magnif, V_prime_ptop_to_pmid*magnif, ...
    U_g_prime_pmid_to_pbot*magnif, V_g_prime_pmid_to_pbot*magnif};

title_chc = {'U''', 'V''', 'U_{g}''', 'V_{g}'''};
z1_chc = {z_top, z_top, z_mid, z_mid};
z2_chc = {z_mid, z_mid, z_bot, z_bot};

for sp = 1 : rowcols(1)*rowcols(2)
    minmax{sp} = [cmaps_cont(1) cmaps_cont(end)];
    cmaps_custom{sp} = cmapcust(cmaps,cmaps_cont);
    
    axis_setup{sp} = [lon_min lon_max lat_min lat_max];
    x{sp} = x_chc{x_ind(sp)};
    y{sp} = y_chc{x_ind(sp)};

end

lett = 'a':'z';
font_size = 13;
fig_color = [0.7 0.7 0.7];

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
set(fig,'units','centimeters',...
    'position',[0 0 ...
    (marg_l+colN*x_sp+gap_w*(colN-1)+marg_r) ...
    (marg_b+rowN*y_sp+gap_h*(rowN-1)+marg_t)]);


for sp = 1 : rowN*colN
    ax = axes('Units','centimeters', ...
        'Position',[...
        (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)), ...
        (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
        x_sp, ...
        y_sp]);
    
    colormap(ax, cmaps_custom{sp});
    pcolor(x{sp}, y{sp}, data{sp})
    axis(axis_setup{sp})
    shading interp
    caxis([minmax{sp}(1) minmax{sp}(2)]);
    hold on
    
    % h = plot(lon_u_ALLC_repelem, lat_v_SBC_centr_repelem, ...
    %     'k--', 'linewidth', 0.5);
    
    if sp <= 2
        %
        lon_ind = lon_u_ALLC_repelem >= ...
            aus8_currents.nudging.north.west_GAB.lon_u_east & ...
            lon_u_ALLC_repelem <= ...
            aus8_currents.nudging.north.eGs_wBS.lon_u_west;
        h = plot(lon_u_ALLC_repelem(lon_ind), ...
            lat_v_SBC_north_repelem(lon_ind), ...
            'g', 'linewidth', 2);
        %
        lon_ind = lon_u_ALLC_repelem >= ...
            aus8_currents.nudging.north.eGs_wBS.lon_u_east & ...
            lon_u_ALLC_repelem <= ...
            aus8_currents.nudging.north.sBS_wT.lon_u_west;
        h = plot(lon_u_ALLC_repelem(lon_ind), ...
            lat_v_SBC_north_repelem(lon_ind), ...
            'g', 'linewidth', 2);
    end
    %
    h = plot(lon_u_ALLC_repelem, lat_v_SBC_south_repelem, ...
        'r', 'linewidth', 2);
    % h = plot(lon_u_ALLC_repelem, lat_v_DRC_centr_repelem, ...
    %     'k-', 'linewidth', 0.5);
    if sp <= 2
        h = plot(lon_u_ALLC_repelem, lat_v_DRC_south_repelem, ...
            '--m', 'linewidth', 2);
    else
        h = plot(lon_u_ALLC_repelem, lat_v_DRC_south_repelem, ...
            ':m', 'linewidth', 2);
    end
    
    
    if sp == 1
        arrow([139 -33], ...
            [lon_u_ALLC_repelem(320) lat_v_SBC_north_repelem(320)], 8)
        text(139, -33, '$y_{CC}$', 'fontsize',font_size)
        arrow([139 -34], ...
            [lon_u_ALLC_repelem(320) lat_v_SBC_south_repelem(320)], 8)
        text(139, -34, '$y_{int}$', 'fontsize',font_size)
        arrow([139 -35], ...
            [lon_u_ALLC_repelem(320) lat_v_DRC_south_repelem(320)], 8)
        text(139, -35, 'upper $y_{OF}$', 'fontsize',font_size)
    end
    
    if sp == 3
        arrow([135 -33.8], ...
            [131.6 -33.8], 8)
        text(135, -33.8, '$z_{int}$ area', 'fontsize',font_size)
        arrow([139 -35], ...
            [lon_u_ALLC_repelem(320) lat_v_DRC_south_repelem(320)], 8)
        text(139, -35, 'lower $y_{OF}$', 'fontsize',font_size)
    end
    
    title(['(' lett(sp) ') $' title_chc{sp} ...
        '$ ($m^{2}/s$) integrated from ' ...
    '$z=' num2str(z1_chc{sp}) '$ to $z=' num2str(z2_chc{sp}) '$ $m$'])
    grid
    set(gca,'layer','top','color',fig_color,...
        'fontsize',font_size,'tickdir','out', ...
        'xtick', lon_min:2:lon_max, 'ytick', lat_min:2:lat_max)
    if row_ind(sp) ~= rowN, set(gca,'xticklabel',''), end
    if col_ind(sp) ~= 1, set(gca,'yticklabel',''), end
    
    if sp == rowN*colN
        ax = axes;
        set(gca,'Visible','off')
        colormap(ax, cmaps);
        cbar = colorbar;
        set(cbar,'ytick',cmaps_linspace, ...
            'YAxisLocation','right','YTickLabel',cmaps_y_label,...
            'fontsize',font_size)
        set(cbar,'units','centimeters','position', [...
            (marg_l+x_sp*(cm(sp))+gap_w*(cm(sp)-1)+plot_cbar_gap), ...
            (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
            cbar_x, ...
            cbar_y]);
    end
end

outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
export_fig(fig, ...
    [figures_path mfilename '/' scriptname(1:3) ...
    '_fig' num2str(fig_n) '_'], ...
    '-transparent', '-m2.5')
close


%%
play(bird_f_path,[1 (get(bird_f_path, 'SampleRate')*3)]);

