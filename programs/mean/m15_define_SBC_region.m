%% define current boundaries
clearvars('-except', '*_path')

load aus8_ZD_method
pres = (0:10:2000)';
pres_mid = aus8_ZD_method.pres_mid;
depth_h_raw = aus8_ZD_method.depth_thicknesses;
lat_u = aus8_ZD_method.lat_u;
lon_u = aus8_ZD_method.lon_u;
u_g_prime = aus8_ZD_method.u_g_prime;
lat_v = aus8_ZD_method.lat_v;
lon_v = aus8_ZD_method.lon_v;
v_g_prime = aus8_ZD_method.v_g_prime;
lat_F = aus8_ZD_method.lat_F;
lon_F = aus8_ZD_method.lon_F;
div_UV_prime = aus8_ZD_method.div_UV_prime;
U_ek = aus8_ZD_method.U_ek;
V_ek = aus8_ZD_method.V_ek;


%% define bottom
%%% BOTTOM PRES
% pressure of the bottom surface currents
SC_bottom_pres = 200;
% pressure above the bottom
SC_bottom_pres_mid = 205;
% index vector of the bottom surface currents
SC_bottom_pres_mid_ind = pres_mid == SC_bottom_pres_mid;
% index vector of currents pressure from surface to interface
SC_bottom_all_pres_mid_ind = pres_mid <= SC_bottom_pres_mid;
%%% BOTTOM PRES


%%% WEST & EAST AND LON UV IND
% longitude of western boundary of the control volume
ALLC_west_lon_u = 115;
% index of longitude of western boundary of the control volume
ALLC_west_lon_u_ind = find(lon_u == ALLC_west_lon_u);
% longitude of eastern boundary of the control volume
ALLC_east_lon_u = 147;
% index of longitude of eastern boundary of the control volume
ALLC_east_lon_u_ind = find(lon_u == ALLC_east_lon_u);
% vector ind
ALLC_west_to_east_lon_u_ind = ALLC_west_lon_u_ind:ALLC_east_lon_u_ind;
aus8_currents.ALLC_west_lon_u = ALLC_west_lon_u;
aus8_currents.ALLC_east_lon_u = ALLC_east_lon_u;
aus8_currents.ALLC_west_to_east_lon_u_ind = ALLC_west_to_east_lon_u_ind;
%%% WEST & EAST LON U


%%% LAT V BOTTOM EXTENT
% southern extent of the surface current box, taken from the first land
% at the interface between the surface and the subsurface currents
SC_v_bottom_south_1 = 0.625; % degree latitude
SC_v_bottom_south_2 = 0.875;
SC_v_bottom_south_3 = 1.125;
SC_v_bottom_south_4 = 1.375;

% length of the southern extent line at the interface between the SC SUBSC
SC_v_bottom_south_length_1 = SC_v_bottom_south_1 / 0.125;
SC_v_bottom_south_length_2 = SC_v_bottom_south_2 / 0.125;
SC_v_bottom_south_length_3 = SC_v_bottom_south_3 / 0.125;
SC_v_bottom_south_length_4 = SC_v_bottom_south_4 / 0.125;
% index vector of the southern extent line at the interface between the SC
% SUBSC
SC_v_bottom_south_vec_1 = 1:SC_v_bottom_south_length_1;
SC_v_bottom_south_vec_2 = 1:SC_v_bottom_south_length_2;
SC_v_bottom_south_vec_3 = 1:SC_v_bottom_south_length_3;
SC_v_bottom_south_vec_4 = 1:SC_v_bottom_south_length_4;
SC_v_bottom_south_length_1 = length(SC_v_bottom_south_vec_1);
SC_v_bottom_south_length_2 = length(SC_v_bottom_south_vec_2);
SC_v_bottom_south_length_3 = length(SC_v_bottom_south_vec_3);
SC_v_bottom_south_length_4 = length(SC_v_bottom_south_vec_4);
%
SC_u_bottom_south_vec_1 = 1:SC_v_bottom_south_length_1-1;
SC_u_bottom_south_vec_2 = 1:SC_v_bottom_south_length_2-1;
SC_u_bottom_south_vec_3 = 1:SC_v_bottom_south_length_3-1;
SC_u_bottom_south_vec_4 = 1:SC_v_bottom_south_length_4-1;
SC_u_bottom_south_length_1 = length(SC_u_bottom_south_vec_1);
SC_u_bottom_south_length_2 = length(SC_u_bottom_south_vec_2);
SC_u_bottom_south_length_3 = length(SC_u_bottom_south_vec_3);
SC_u_bottom_south_length_4 = length(SC_u_bottom_south_vec_4);
%%% LAT V BOTTOM EXTENT


%%% WEST TO EAST LON U AND LON V
% place the u longitudes within the west and east boundaries into the
% template
lon_u_SC_bottom_line = lon_u(ALLC_west_to_east_lon_u_ind);
% lon u grid
SC_bottom_lon_u = ...
    repmat(lon_u_SC_bottom_line, [SC_u_bottom_south_length_4, 1]);
% lon v
lon_v_SC_bottom_line = lon_v(ALLC_west_to_east_lon_u_ind(1:end-1));
SC_bottom_lon_v = ...
    repmat(lon_v_SC_bottom_line,[SC_v_bottom_south_length_4, 1]);
%%% WEST TO EAST LON U AND LON V


%%% PREP FOR LAT V AT THE BOTTOM
% templates for:
% latitudes of the interface lines
% indices vectors of the latitudes of the interface lines
SC_bottom_lat_u = NaN(...
    SC_u_bottom_south_length_4, ...
    length(lon_u_SC_bottom_line));

SC_bottom_lat_v = NaN(...
    SC_v_bottom_south_length_4, ...
    length(lon_v_SC_bottom_line));

SC_bottom_lat_v_south_line = ...
    NaN(1, length(lon_v_SC_bottom_line));

% longitude index counter
jj_count = 0;
%%% PREP FOR LAT V AT THE BOTTOM


%%% GET LAT U AND V AT THE BOTTOM
% for all longitudes within the west and east boundaries
for jj = ALLC_west_to_east_lon_u_ind
    % increase lon counter
    jj_count = jj_count + 1;
    
    % lat v
    SC_bottom_all_lat_u = u_g_prime(:,jj,SC_bottom_pres_mid_ind);
    SC_bottom_lat_u_first_ocean_ind = ...
        find(isnan(SC_bottom_all_lat_u), 1, 'last')+1;
    
    if lon_u(jj) <= 134
        SC_u_bottom_south_vec = SC_u_bottom_south_vec_1;
        SC_v_bottom_south_vec = SC_v_bottom_south_vec_1;
    elseif lon_u(jj) > 134 && lon_u(jj) <= 142
        SC_u_bottom_south_vec = SC_u_bottom_south_vec_2;
        SC_v_bottom_south_vec = SC_v_bottom_south_vec_2;
    elseif lon_u(jj) > 142 && lon_u(jj) <= 142.25
        SC_u_bottom_south_vec = SC_u_bottom_south_vec_3;
        SC_v_bottom_south_vec = SC_v_bottom_south_vec_3;
    elseif lon_u(jj) > 142.25 && lon_u(jj) <= 143.25
        SC_u_bottom_south_vec = SC_u_bottom_south_vec_4;
        SC_v_bottom_south_vec = SC_v_bottom_south_vec_4;
    elseif lon_u(jj) > 143.25 && lon_u(jj) <= 144.75
        SC_u_bottom_south_vec = SC_u_bottom_south_vec_3;
        SC_v_bottom_south_vec = SC_v_bottom_south_vec_3;
    elseif lon_u(jj) > 144.75 && lon_u(jj) <= 145.875
        SC_u_bottom_south_vec = SC_u_bottom_south_vec_2;
        SC_v_bottom_south_vec = SC_v_bottom_south_vec_2;
    elseif lon_u(jj) > 145.875
        SC_u_bottom_south_vec = SC_u_bottom_south_vec_1;
        SC_v_bottom_south_vec = SC_v_bottom_south_vec_1;
    end
    
    SC_bottom_lat_u_ind_now = ...
        (SC_bottom_lat_u_first_ocean_ind+SC_u_bottom_south_vec-1)';
    SC_bottom_lat_u(SC_u_bottom_south_vec,...
        jj_count) = lat_u(SC_bottom_lat_u_ind_now);
    
    if jj ~= ALLC_west_to_east_lon_u_ind(end)
        SC_bottom_lat_v_ind_now = ...
            (SC_bottom_lat_u_first_ocean_ind+SC_v_bottom_south_vec-1)';
        %
        SC_bottom_lat_v(SC_v_bottom_south_vec,...
            jj_count) = lat_v(SC_bottom_lat_v_ind_now);
        
        SC_bottom_lat_v_south_line(jj_count) = ...
            lat_v(SC_bottom_lat_v_ind_now(end));
    end
    
end
%%% GET LAT U AND V AT THE BOTTOM


%% Define surface
%
SC_u = zeros(size(u_g_prime));
SC_v = zeros(size(v_g_prime));
SC_pres_mid_length = length(find(pres_mid < SC_bottom_pres));

%
jj_count = 0;
for jj = ALLC_west_to_east_lon_u_ind
    jj_count = jj_count +1;
    
    %
    SC_bottom_lat_u_ind_now = ...
        find(ismember(lat_u, SC_bottom_lat_u(:,jj_count)));
    %
    SC_u(...
        SC_bottom_lat_u_ind_now, ...
        jj, ...
        SC_pres_mid_length(end)) = u_g_prime(...
        SC_bottom_lat_u_ind_now, ...
        jj, ...
        SC_pres_mid_length(end));
    
    if jj ~= ALLC_west_to_east_lon_u_ind(end)
        %
        SC_bottom_lat_v_ind_now = ...
            find(ismember(lat_v, SC_bottom_lat_v(:,jj_count)));
        %
        SC_v(...
            SC_bottom_lat_v_ind_now, ...
            jj, ...
            SC_pres_mid_length(end)) = v_g_prime(...
            SC_bottom_lat_v_ind_now, ...
            jj, ...
            SC_pres_mid_length(end));
    end
end


%%% GET U AND V
for kk = SC_pres_mid_length-1:-1:1
    % longitude index counter
    jj_count = 0;
    % for all longitudes within the west and east boundaries
    for jj = ALLC_west_to_east_lon_u_ind
        % increase lon counter
        jj_count = jj_count + 1;
        
        % lat u
        SC_lat_u_now = u_g_prime(:,jj,kk);
        
        if jj >= 57 && jj <= 65
            SC_lat_u_now(1:32) = NaN;
        end
        
        
        flag = 0;
        jj_1 = 129;
        if jj >= 57 && jj <= jj_1
            p_now = 15; p_now_ind = find(pres_mid == p_now)+1;
            if pres_mid(kk) <= p_now
                flag = 1;
                mask_ind = SC_v(:,jj,p_now_ind) == 0;
                SC_lat_u_now(mask_ind) = NaN;
            end
            
            jj_2 = 138;
        elseif jj > jj_1 && jj <= jj_2
            p_now = 105; p_now_ind = find(pres_mid == p_now)+1;
            if pres_mid(kk) <= p_now
                flag = 1;
                mask_ind = SC_v(:,jj,p_now_ind) == 0;
                SC_lat_u_now(mask_ind) = NaN;
            end
            
        elseif jj > jj_2
            p_now = 105; p_now_ind = find(pres_mid == p_now)+1;
            if pres_mid(kk) <= p_now
                flag = 1;
                mask_ind = SC_v(:,jj,p_now_ind) == 0;
                SC_lat_u_now(mask_ind) = NaN;
            end
        end
        
        SC_lat_u_first_ocean_now_ind = ...
            find(isfinite(SC_lat_u_now), 1, 'first');
        
        
        if jj ~= ALLC_west_to_east_lon_u_ind(end)
            % south
            SC_lat_u_south_ind_now = ...
                find(ismember(lat_u, ...
                SC_bottom_lat_v_south_line(jj_count)+1/16));
            
            SC_lat_u_north_south_vec_ind_now = ...
                (SC_lat_u_first_ocean_now_ind : ...
                SC_lat_u_south_ind_now)';
            
            %
            SC_lat_v_north_south_vec_ind_now = ...
                SC_lat_u_north_south_vec_ind_now;
            SC_lat_v_north_south_vec_ind_now(end+1) = ...
                SC_lat_u_north_south_vec_ind_now(end)+1; %#ok<SAGROW>
            v_g_prime_now = v_g_prime(...
                SC_lat_v_north_south_vec_ind_now, ...
                jj, ...
                kk);
            SC_v(...
                SC_lat_v_north_south_vec_ind_now, ...
                jj, ...
                kk) = v_g_prime_now;
        end
        
        %
        u_g_prime_now = u_g_prime(...
            SC_lat_u_north_south_vec_ind_now, ...
            [jj, jj+1], ...
            kk);
        SC_u(...
            SC_lat_u_north_south_vec_ind_now, ...
            [jj, jj+1], ...
            kk) = u_g_prime_now;
        
    end
end
%%% GET U AND V

SC_u_ind = SC_u ~= 0;
SC_v_ind = SC_v ~= 0;


% GET NORTH SOUTH BOUNDARIES FOR CONTOUR INTEGRATION
% lon u
SC_lon_u_all = repmat(lon_u, [length(lat_u), 1, length(pres_mid)]);
SC_lon_u_all(~SC_u_ind) = NaN;
SC_lon_u = repmat(...
    lon_u(ALLC_west_to_east_lon_u_ind), ...
    [length(pres_mid), 1]);

% lat v
SC_lat_v_all = repmat(lat_v, [1, length(lon_v), length(pres_mid)]);
SC_lat_v_all(~SC_v_ind) = NaN;
[SC_lat_v_north, SC_lat_v_south] = ...
    deal(NaN(length(pres_mid), length(ALLC_west_to_east_lon_u_ind)-1));
for kk = 1 : length(pres_mid)
    jj_count = 0;
    for jj = ALLC_west_to_east_lon_u_ind
        jj_count = jj_count + 1;
        SC_lat_v_north_ind_now = ...
            find(isfinite(SC_lat_v_all(:,jj,kk)), 1, 'first');
        SC_lat_u_south_ind_now = ...
            find(isfinite(SC_lat_v_all(:,jj,kk)), 1, 'last');
        
        if ~isempty(SC_lat_v_north_ind_now)
            SC_lat_v_north(kk, jj_count) = ...
                SC_lat_v_all(SC_lat_v_north_ind_now, jj, kk);
            SC_lat_v_south(kk, jj_count) = ...
                SC_lat_v_all(SC_lat_u_south_ind_now, jj, kk);
        end
    end
end

%
SC_lon_u_repelem = [...
    SC_lon_u(:,1), ...
    repelem(SC_lon_u(:,2:end-1), 1, 2), ...
    SC_lon_u(:,end)];
SC_lat_v_north_repelem = ...
    repelem(SC_lat_v_north, 1, 2);
SC_lat_v_south_repelem = ...
    repelem(SC_lat_v_south, 1, 2);


%%% SAVE
aus8_currents.SC.u_g_prime_ind = SC_u_ind;
aus8_currents.SC.v_g_prime_ind = SC_v_ind;
aus8_currents.SC.lon_u = SC_lon_u;
aus8_currents.SC.lat_v_north = SC_lat_v_north;
aus8_currents.SC.lat_v_south = SC_lat_v_south;
aus8_currents.SC.lon_u_repelem = SC_lon_u_repelem;
aus8_currents.SC.lat_v_north_repelem = SC_lat_v_north_repelem;
aus8_currents.SC.lat_v_south_repelem = SC_lat_v_south_repelem;
%%% SAVE


%% plot prep
%
depth_h = permute(depth_h_raw, [3 2 1]);
depth_h_u = repmat(depth_h, [length(lat_u), length(lon_u)]);
depth_h_v = repmat(depth_h, [length(lat_v), length(lon_v)]);

pres_mid_now = 5;
if pres_mid_now == 5
    quiv_S = 10;
else
    quiv_S = 3;
end

p_ind = find(pres_mid == pres_mid_now);

%
lon_plot_west = 114;
lon_plot_east = 148;
lon_F_ind = find(lon_F >= lon_plot_west & lon_F <= lon_plot_east);
lon_u_ind = find(lon_u >= lon_plot_west & lon_u <= lon_plot_east);
lon_v_ind = find(lon_v >= lon_plot_west & lon_v <= lon_plot_east);
lat_plot_north = -33;
lat_plot_south = -45;
lat_F_ind = find(lat_F <= lat_plot_north & lat_F >= lat_plot_south);
lat_u_ind = find(lat_u <= lat_plot_north & lat_u >= lat_plot_south);
lat_v_ind = find(lat_v <= lat_plot_north & lat_v >= lat_plot_south);

%
data.SC_u = SC_u(lat_u_ind,lon_u_ind,p_ind);
data.SC_v = SC_v(lat_v_ind,lon_v_ind,p_ind);
data.div_UV_prime = div_UV_prime;


%% plot maps of U and V SBC
close all
fig1 = figure(1);

set(gcf,'units','normalized','outerposition',[0 0 0.9 0.9])

rowN = 1; colN = 1;
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (flipud(repmat((1:rowN)',1,colN)))';
gap_w = 0.02; % gap width between subplots
gap_h = 0.06; % gap height between subplots
marg_b = 0.03; % bottom margin
marg_t = 0.03; % top margin
marg_l = 0.02; % left margin
marg_r = 0.005; % right margin
x_sp = (1 - marg_l - marg_r - gap_w*(colN-1))/colN; % x subplot length
y_sp = (1 - marg_b - marg_t - gap_h*(rowN-1))/rowN; % y subplot length

h_axes_sp = tight_subplot(rowN, colN, ...
    [gap_h gap_w], [marg_b marg_t], [marg_l marg_r]);

axes(h_axes_sp(1))

colormap([1 1 1]);
h_pcolor = pcolor(...
    lon_F(lon_F_ind)-0.0625, ...
    lat_F(lat_F_ind)+0.0625, ...
    data.div_UV_prime(lat_F_ind,lon_F_ind));


hold on
% west
plot([SC_lon_u(1,1) SC_lon_u(1,1)], ...
    [SC_lat_v_north(p_ind,1) SC_lat_v_south(p_ind,1)], ...
    'g','linewidth',0.5,'linestyle','-')
% north
plot(SC_lon_u_repelem(p_ind,:), ...
    SC_lat_v_north_repelem(p_ind,:), ...
    'g','linewidth',0.5,'linestyle','-')
% south
plot(SC_lon_u_repelem(p_ind,:), ...
    SC_lat_v_south_repelem(p_ind,:), ...
    'g','linewidth',0.5,'linestyle','-')
% east
plot([SC_lon_u(1,end) SC_lon_u(1,end)], ...
    [SC_lat_v_north(p_ind,end) SC_lat_v_south(p_ind,end)], ...
    'g','linewidth',0.5,'linestyle','-')

hold on
[lat_u_g_prime_NaN_ind, lon_u_g_prime_NaN_ind] = ...
    find(isnan(u_g_prime(lat_u_ind,lon_u_ind,p_ind)));

lat_u_zoom = lat_u(lat_u_ind);
lon_u_zoom = lon_u(lon_u_ind);

lat_u_g_prime_NaN = lat_u_zoom(lat_u_g_prime_NaN_ind);
lon_u_g_prime_NaN = lon_u_zoom(lon_u_g_prime_NaN_ind);

scatter(lon_u_g_prime_NaN, lat_u_g_prime_NaN, 4, ...
    'o', 'k');

hold on
[lat_v_g_prime_NaN_ind, lon_v_g_prime_NaN_ind] = ...
    find(isnan(v_g_prime(lat_v_ind,lon_v_ind,p_ind)));

lat_v_zoom = lat_v(lat_v_ind);
lon_v_zoom = lon_v(lon_v_ind);

lat_v_g_prime_NaN = lat_v_zoom(lat_v_g_prime_NaN_ind);
lon_v_g_prime_NaN = lon_v_zoom(lon_v_g_prime_NaN_ind);

scatter(lon_v_g_prime_NaN, lat_v_g_prime_NaN, 4, ...
    'o', 'k');

hold on
nn = 1;
quiv_S = 5;
SBC_U_g_times_depth = SC_u.*depth_h_u;
SBC_U_g = nansum(SBC_U_g_times_depth,3);
SBC_U_g_mask = SBC_U_g==0;
SBC_U0 = SBC_U_g + U_ek;
SBC_U0(SBC_U_g_mask) = 0;
SBC_U = SBC_U0(lat_u_ind,lon_u_ind);
[lon_mg, lat_mg]=meshgrid(lon_u(lon_u_ind),lat_u(lat_u_ind));
v_g_0 = zeros(size(SBC_U));
quiver(...
    lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
    SBC_U(1:nn:end, 1:nn:end), ...
    v_g_0(1:nn:end, 1:nn:end), ...
    quiv_S, 'r');

SBC_V_g_times_depth = SC_v.*depth_h_v;
SBC_V_g_raw = nansum(SBC_V_g_times_depth,3);
SBC_V_g_raw_mask = SBC_V_g_raw==0;
SBC_V_raw0 = SBC_V_g_raw + V_ek;
SBC_V_raw0(SBC_V_g_raw_mask) = 0;
SBC_V_raw = SBC_V_raw0(lat_v_ind,lon_v_ind);

SBC_V = NaN(size(SBC_V_raw));
for jj = 1 : length(lon_v_zoom)
    first_now = find(SBC_V_raw(:,jj)~=0,1,'first');
    last_now = find(SBC_V_raw(:,jj)~=0,1,'last');
    SBC_V([first_now,last_now],jj) = SBC_V_raw([first_now,last_now],jj);
end

[lon_mg, lat_mg]=meshgrid(lon_v(lon_v_ind),lat_v(lat_v_ind));
u_g_0 = zeros(size(SBC_V));
quiver(...
    lon_mg(1:nn:end, 1:nn:end), lat_mg(1:nn:end, 1:nn:end), ...
    u_g_0(1:nn:end, 1:nn:end), ...
    SBC_V(1:nn:end, 1:nn:end), ...
    quiv_S, 'b');


title(['Fluxes at ' num2str(pres_mid_now) ' dbar'])


grid
font_size = 10;
set(gca,'layer','top','color',[0.7 0.7 0.7],...
    'fontsize',font_size,'fontweight','bold')

outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
export_fig(fig1, [figures_path mfilename '/' scriptname(1:3) '_'], ...
    '-m4')
close


%% save
save([cars_out_path 'aus8_ZD_method'], 'aus8_ZD_method')
disp(['aus8_ZD_method saved in ' ...
    cars_out_path 'aus8_ZD_method'])


