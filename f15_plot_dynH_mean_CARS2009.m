%%
clearvars('-except', '*_path')
play(bird_i_path,[1 (get(bird_i_path, 'SampleRate')*3)]);

load([data_path 'SACS_data/hgt2000_cars2009'])
% temperature
[temp_data, temp_gen_att, temp_att] = ...
    nc2mat([data_path ...
    'KDS75_cp/seasonal_climatology.nc'], ...
    'ALL');

%%
days_in_months = [...
    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
time_step = NaN(1,12);
time_step(1) = days_in_months(1)/2;
for p = 2 : 12
    time_step(p) = ...
        time_step(p-1) + days_in_months(p-1)/2 + days_in_months(p)/2;
end

Months = {...
    'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

% time cycle setup
time_serie = 2 * pi / 365 * time_step;

for ll = 1 : 12
    temp_monthly = ...
        ... % mean component
        Var.mean + ...
        ... % annual component
        Var.an_cos*cos(time_serie(ll)) + ...
        Var.an_sin*sin(time_serie(ll)) + ...
        ... % semi-annual component
        Var.sa_cos*cos(2*time_serie(ll)) + ...
        Var.sa_sin*sin(2*time_serie(ll));
    
    aus8_monthly.(Months{ll}) = temp_monthly;
    disp([Months{ll} ' OK!'])
end


JFM = aus8_monthly.Jan;
JFM(:,:,2) = aus8_monthly.Feb;
JFM(:,:,3) = aus8_monthly.Mar;
hgt2000 = mean(JFM,3)';

JFM = aus8_monthly.Apr;
JFM(:,:,2) = aus8_monthly.May;
JFM(:,:,3) = aus8_monthly.Jun;
hgt2000(:,:,2) = mean(JFM,3)';

JFM = aus8_monthly.Jul;
JFM(:,:,2) = aus8_monthly.Aug;
JFM(:,:,3) = aus8_monthly.Sep;
hgt2000(:,:,3) = mean(JFM,3)';

JFM = aus8_monthly.Oct;
JFM(:,:,2) = aus8_monthly.Nov;
JFM(:,:,3) = aus8_monthly.Dec;
hgt2000(:,:,4) = mean(JFM,3)';

hgt2000 = (Var.mean - nanmean(nanmean(Var.mean))*1.18)';

sea_level = mean(permute(temp_data.sea_level, ([2 1 3])), 3);


%% 5) plot maps of U and V SBC
Seasons = {'summer', 'summer', 'autumn', 'autumn', ...
    'winter', 'winter','spring', 'spring'};

screen_ratio = 0.75;
fig_n = 1;
rowcols = [2 1];
rowcols_size = [14 6]/screen_ratio/2; % cm
margs = [0.9 0.2 1.8 0.6]/screen_ratio; % cm
gaps = [0.4 0.8]/screen_ratio; % cm
plot_cbar_gap = 1/screen_ratio;
cbar_x = rowcols_size(1);
cbar_y = 0.2/screen_ratio;

% lon_min = 115; lon_max = 147; lat_min = -47; lat_max = -32;
lon_min = 110; lon_max = 152; lat_min = -48; lat_max = -31;

x_chc = {Var.lon, temp_data.xt_ocean+360};
x_ind = [1 2];
y_chc = {Var.lat, temp_data.yt_ocean};

magnif = 1000;
cmap1_cont = -[0.3 0.2 0.15 0.1 0.075 0.05 0.025 0.01 0.001 0]*magnif;
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
cc = 0;
for t = 1
    cc = cc + 1;
    data{t} = hgt2000(:,:,cc)*magnif;
    minmax{t} = [cmaps_cont(1) cmaps_cont(end)];
    cmaps_custom{t} = cmapcust(cmaps,cmaps_cont);
    axis_setup{t} = [lon_min lon_max lat_min lat_max];
    x{t} = x_chc{x_ind(1)};
    y{t} = y_chc{x_ind(1)};
    U_title{t} = 'CARS2009';
end

cc = 0;
for t = 2
    cc = cc + 1;
    data{t} = sea_level(:,:,cc)*magnif;
    minmax{t} = [cmaps_cont(1) cmaps_cont(end)];
    cmaps_custom{t} = cmapcust(cmaps,cmaps_cont);
    axis_setup{t} = [lon_min lon_max lat_min lat_max];
    x{t} = x_chc{x_ind(2)};
    y{t} = y_chc{x_ind(2)};
    U_title{t} = 'MOM01-75z';

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
    
    colormap(ax, cmaps_custom{sp});
    pcolor(x{sp}, y{sp}, data{sp})
    axis(axis_setup{sp})
    shading interp
    caxis([minmax{sp}(1) minmax{sp}(2)]);
    
%     title_chc = ...
%     {'{\boldmath{$V_{up}$}} (arrows) and {$U_{up}$} (shadings)', ...
%     '{\boldmath{$V_{low}$}} (arrows) and {$U_{low}$} (shadings)'};
    
    h_tit = title(['(' lett(sp) ') ' U_title{sp} ' ' Seasons{sp} ...
        ' sea level anomaly'], ...
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
    if col_ind(sp) ~= 1
        set(gca,'yticklabel','')
    else
        ylabel('Latitude')
    end
%     
%     if ~mod(sp,2)
    hold on
    plot([110, 152], [-42, -42], 'k:')
%     end
    
    if sp == rowN*colN
        ax = axes('visible', 'off');
        colormap(ax, cmaps);
        cbar = colorbar('horizontal');
        set(cbar,'ytick',cmaps_linspace, ...
            'YAxisLocation','right','YTickLabel',cmaps_y_label*magnif,...
            'fontsize',font_size,'ticklength',cbar_tick_length)
        set(cbar,'units','centimeters','position', [...
            (marg_l+x_sp*(cm(sp)-2)+gap_w*(cm(sp)-2)), ...
            (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)-plot_cbar_gap), ...
            cbar_x*2+gap_w, ...
            cbar_y]);
        set(get(cbar,'xlabel'),'String','$\eta$ ($mm$)', ...
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
