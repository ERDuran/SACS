%% calculate sc transport
clearvars('-except', 'outputpath')
load aus8_ZD_method
load aus8_currents


%%
% lat lon pres
lat_u = aus8_ZD_method.lat_u;
lon_u = aus8_ZD_method.lon_u;
lat_v = aus8_ZD_method.lat_v;
lon_v = aus8_ZD_method.lon_v;
pres = aus8_ZD_method.pres_mid;

% for dx dy calcs
a = aus8_ZD_method.a.value;
pi180 = aus8_ZD_method.a.pi180;
f = aus8_ZD_method.f.cst_lat;

% thickness for dz calcs
depth_h_raw = aus8_ZD_method.depth_thicknesses;
depth_h = permute(depth_h_raw, [3 2 1]);
depth_h_u = repmat(depth_h, [length(lat_u), length(lon_u)]);
depth_h_v = repmat(depth_h, [length(lat_v), length(lon_v)]);

% dx for v and dy for u calcs
dx_v = NaN(length(lat_v), length(lon_u)-1);
for ii = 1 : length(lat_v)
    dx_v(ii,:) = a * cos(lat_v(ii) * pi180) * ...
        (lon_u(2:end) - lon_u(1:end-1)) * pi180;
end
dy_raw = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
dy_u = repmat(dy_raw, [1 length(lon_u)]);


% original uv g prime
u_g_prime = aus8_ZD_method.u_g_prime;
v_g_prime = aus8_ZD_method.v_g_prime;
u_g_prime(isnan(u_g_prime)) = 0;
v_g_prime(isnan(v_g_prime)) = 0;

% index of SC uv g prime
SC_u_g_prime_ind = aus8_currents.SC.u_g_prime_ind;
SC_v_g_prime_ind = aus8_currents.SC.v_g_prime_ind;

% SC uv g prime
SC_u_g_prime_all = NaN(size(u_g_prime));
SC_v_g_prime_all = NaN(size(v_g_prime));
SC_u_g_prime_all(SC_u_g_prime_ind) = u_g_prime(SC_u_g_prime_ind);
SC_v_g_prime_all(SC_v_g_prime_ind) = v_g_prime(SC_v_g_prime_ind);

% repelems
SC_lon_u = aus8_currents.SC.lon_u(1,:);
SC_lon_v = SC_lon_u(1:end-1) + 1/16; 
SC_lat_v_north = aus8_currents.SC.lat_v_north(1,:);
SC_lat_v_south = aus8_currents.SC.lat_v_south(1,:);
SC_lon_u_repelem = aus8_currents.SC.lon_u_repelem(1,:);
SC_lat_v_north_repelem = aus8_currents.SC.lat_v_north_repelem(1,:);
SC_lat_v_south_repelem = aus8_currents.SC.lat_v_south_repelem(1,:);

% UV ekman
U_ek = aus8_ZD_method.U_ek;
V_ek = aus8_ZD_method.V_ek;
U_ek(isnan(U_ek)) = 0;
V_ek(isnan(V_ek)) = 0;

% div(U',V')
div_UV_prime = aus8_ZD_method.div_UV_prime;


% lon index of domain
ALLC_lon_u_ind = aus8_currents.ALLC_west_to_east_lon_u_ind;


%% U prime trans
% U prime transport (non-cumulated)
SC_u_g_prime_times_depth_h_u = SC_u_g_prime_all .* depth_h_u;
SC_U_g_prime = nansum(SC_u_g_prime_times_depth_h_u, 3);
SC_u_g_prime_surf_isnan_ind = isnan(SC_u_g_prime_all(:,:,1));
SC_U_g_prime(SC_u_g_prime_surf_isnan_ind) = NaN;
SC_U_prime = SC_U_g_prime + U_ek;
SC_U_trans_mcps_indiv = SC_U_prime .* dy_u;


%% U within
SC_U_trans_mcps_indiv_within = NaN(size(SC_U_prime));

% SC_U_trans_mcps_indiv_within(37:44,57) = ...
%     SC_U_trans_mcps_indiv(37:44,57);
% SC_U_trans_mcps_indiv_within(113:117,313) = ...
%     SC_U_trans_mcps_indiv(113:117,313);

ALLC_lon_u_ind_NB = ALLC_lon_u_ind(2:end-1);

jj_count = 0;
for jj = 2 : 2 : length(SC_lon_u_repelem(1:end-1))
    jj_count = jj_count + 1;
    
    lat_north_ind_1 = ...
        find(ismember(lat_v, SC_lat_v_north_repelem(jj)));
    lat_north_ind_2 = ...
        find(ismember(lat_v, SC_lat_v_north_repelem(jj+1)));
    lat_north_12 = [lat_north_ind_1, lat_north_ind_2];
    
    lat_south_ind_1 = ...
        find(ismember(lat_v, SC_lat_v_south_repelem(jj)));    
    lat_south_ind_2 = ...
        find(ismember(lat_v, SC_lat_v_south_repelem(jj+1)));
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
        lat_north_12(north_12_ind):lat_south_12(south_12_ind)-1;    
    
    SC_U_trans_mcps_indiv_within(lat_vec_now, ...
        ALLC_lon_u_ind_NB(jj_count)) = ...
        SC_U_trans_mcps_indiv(lat_vec_now,...
        ALLC_lon_u_ind_NB(jj_count));

end

% U SBC
SC_U_trans_mcps = nansum(SC_U_trans_mcps_indiv_within, 1);
SC_U_trans = SC_U_trans_mcps .* 10^-6;


%% V SA AND V SBC
% V prime transport (non-cumulated and cumulated)
SC_v_g_prime_times_depth_h_v = SC_v_g_prime_all .* depth_h_v;
SC_V_g_prime = nansum(SC_v_g_prime_times_depth_h_v, 3);
SC_v_g_prime_surf_isnan_ind = isnan(SC_v_g_prime_all(:,:,1));
SC_V_g_prime(SC_v_g_prime_surf_isnan_ind) = NaN;
SC_V_prime = SC_V_g_prime + V_ek;
SC_V_trans_mcps_indiv = SC_V_prime .* dx_v;


%% get u lines inside
SC_V_trans_mcps_indiv_within = NaN(size(SC_V_prime));

SC_lon_v_ind = find(ismember(lon_v, SC_lon_v));

jj_count = 0;
for jj = SC_lon_v_ind
    jj_count = jj_count + 1;
    lat_v_within_ind = ...
        find(ismember(lat_v, ...
        SC_lat_v_north(jj_count)-1/8:-1/8:SC_lat_v_south(jj_count)+1/8)); 
    SC_V_trans_mcps_indiv_within(lat_v_within_ind, jj) = ...
        SC_V_trans_mcps_indiv(lat_v_within_ind,jj);
end

% U SBC
SC_V_trans_mcps = nansum(SC_V_trans_mcps_indiv_within, 1);
SC_V_trans = SC_V_trans_mcps .* 10^-6;


%% V SC
% V SA
SC_north_trans_nc_mcps = NaN(1, length(lon_v));
for jj = 1 : length(SC_lon_u_repelem(1:end-1))
    lon_u_repelem_ind = find(lon_u == SC_lon_u_repelem(jj));
    
    if mod(jj,2) % if odd
        % then this is a v point, always one and sign doesn't change
        % for V_SA northward is outward, hence southward changes 
        % to positive
        
        SC_V_north_trans_lat_ind = find(isfinite(...
            SC_V_trans_mcps_indiv(:,lon_u_repelem_ind)), 1, 'first');
        
        SC_V_north_trans_now = ...
            SC_V_trans_mcps_indiv(...
            SC_V_north_trans_lat_ind, lon_u_repelem_ind);
        
        SC_V_north_trans = -SC_V_north_trans_now;
       
    else % even
        % then this is a u point, there can be zero or many and the sign
        % can change. The sign depens on the direction of the boundary
        % going west to east, if the boundary is going north then U_SA 
        % eastward is inward, hence eastward stays positive.
        % if the boundary is going southward, then U_SA eastward is outward
        % hence it changes to negative
        
        % but first, if the boundary is not going north or south, then 
        % there is no zonal flow
        
        
        if SC_lat_v_north_repelem(jj) == SC_lat_v_north_repelem(jj+1)
            SC_U_north_trans = 0;
            
        else
            % if boundary is going northward:
            if SC_lat_v_north_repelem(jj) < SC_lat_v_north_repelem(jj+1)
                lat_u_within_ind = find(...
                    lat_u > SC_lat_v_north_repelem(jj) & ...
                    lat_u < SC_lat_v_north_repelem(jj+1));
                
                SC_U_north_trans_now = ...
                    SC_U_trans_mcps_indiv(lat_u_within_ind, lon_u_repelem_ind);
                SC_U_north_trans_all = SC_U_north_trans_now;
                
                
            else % if boundary is going southward
                lat_u_within_ind = find(...
                    lat_u > SC_lat_v_north_repelem(jj+1) & ...
                    lat_u < SC_lat_v_north_repelem(jj));
                
                SC_U_north_trans_now = ...
                    SC_U_trans_mcps_indiv(lat_u_within_ind, lon_u_repelem_ind);
                SC_U_north_trans_all = -SC_U_north_trans_now;                
                
            end
            
            % sum u trans up in case there are more than 1
            SC_U_north_trans = sum(SC_U_north_trans_all);
            
        end
    end
    
    if jj == 1
        SC_north_trans_nc_mcps(57) = SC_V_north_trans;
    elseif mod(jj,2)
        SC_north_trans_nc_mcps(lon_u_repelem_ind) = ...
            SC_V_north_trans + SC_U_north_trans;
    end
end


% V SBC
SC_south_trans_nc_mcps = NaN(1, length(lon_u));
for jj = 1 : length(SC_lon_u_repelem(1:end-1))
    lon_u_repelem_ind = find(lon_u == SC_lon_u_repelem(jj));
    
    if mod(jj,2) % if odd
        % then this is a v point, always one and sign doesn't change
        % for V_SA northward is inward, hence southward stays negative
        
        SC_V_south_trans_lat_ind = find(isfinite(...
            SC_V_trans_mcps_indiv(:,lon_u_repelem_ind)), 1, 'last');
        
        SC_V_south_trans_now = ...
            SC_V_trans_mcps_indiv(...
            SC_V_south_trans_lat_ind, lon_u_repelem_ind);
        
        SC_V_south_trans = SC_V_south_trans_now;
       
    else % even
        % then this is a u point, there can be zero or many and the sign
        % can change. The sign depens on the direction of the boundary
        % going west to east, if the boundary is going north then U_SA 
        % eastward is outward, hence eastward changes to negative.
        % if the boundary is going southward, then U_SA eastward is inward
        % hence it stays positive
        
        % but first, if the boundary is not going north or south, then 
        % there is no zonal flow
        
        
        if SC_lat_v_south_repelem(jj) == SC_lat_v_south_repelem(jj+1)
            SC_U_south_trans = 0;
            
        else
            % if boundary is going northward:
            if SC_lat_v_south_repelem(jj) < SC_lat_v_south_repelem(jj+1)
                lat_u_within_ind = find(...
                    lat_u > SC_lat_v_south_repelem(jj) & ...
                    lat_u < SC_lat_v_south_repelem(jj+1));
                
                SC_U_south_trans_now = ...
                    SC_U_trans_mcps_indiv(lat_u_within_ind, lon_u_repelem_ind);
                SC_U_south_trans_all = -SC_U_south_trans_now;
                
                
            else % if boundary is going southward
                lat_u_within_ind = find(...
                    lat_u > SC_lat_v_south_repelem(jj+1) & ...
                    lat_u < SC_lat_v_south_repelem(jj));
                
                SC_U_south_trans_now = ...
                    SC_U_trans_mcps_indiv(lat_u_within_ind, lon_u_repelem_ind);
                SC_U_south_trans_all = SC_U_south_trans_now;                
                
            end
            
            % sum u trans up in case there are more than 1
            SC_U_south_trans = sum(SC_U_south_trans_all);
            
        end
    end
    
    if jj == 1
        SC_south_trans_nc_mcps(57) = SC_V_south_trans;
    elseif mod(jj,2)
        SC_south_trans_nc_mcps(lon_u_repelem_ind) = ...
            SC_V_south_trans + SC_U_south_trans;
    end
end


% V SBC
SC_north_trans_nc = SC_north_trans_nc_mcps .* 10^-6;
SC_south_trans_nc = SC_south_trans_nc_mcps .* 10^-6;

[SC_north_trans, SC_south_trans] = deal(NaN(1, length(lon_v)));

SC_north_trans(57) = SC_north_trans_nc(57);
SC_south_trans(57) = SC_south_trans_nc(57);

for jj = ALLC_lon_u_ind(2:end)
    SC_north_trans(jj) = ...
        SC_north_trans(jj-1) + SC_north_trans_nc(jj);
    
    SC_south_trans(jj) = ...
        SC_south_trans(jj-1) + SC_south_trans_nc(jj);
    
end


%% SBC div UV prime
%
dx_u = NaN(length(lat_u), length(lon_u)-1);
for ii = 1 : length(lat_u)
    dx_u(ii,:) = a * cos(lat_u(ii) * pi180) * ...
        (lon_u(2:end) - lon_u(1:end-1)) * pi180;
end
dy_raw = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
dy_v = repmat(dy_raw, [1 length(lon_v)]);

du = SC_U_prime(:,2:end) - SC_U_prime(:,1:end-1);
dudx = du ./ dx_u;

lat_v_repmat = repmat(lat_v,1,length(lon_v));
dv = SC_V_prime(1:end-1,:).*cos(lat_v_repmat(1:end-1,:) * pi180) - ...
    SC_V_prime(2:end,:).*cos(lat_v_repmat(2:end,:) * pi180);
lat_u_repmat = repmat(lat_u,1,length(lon_v));
dvdy = 1./cos(lat_u_repmat * pi180).*dv./dy_v;

SC_div_UV_prime = dudx + dvdy;

for ii = 1 : length(lat_u)
    for jj = 1 : length(lon_v)
        if isfinite(SC_div_UV_prime(ii,jj))
            if SC_div_UV_prime(ii,jj) == div_UV_prime(ii,jj)
                SC_div_UV_prime(ii,jj) = 0;
            end
        end
    end
end


%% W SBC
SC_w_star = SC_div_UV_prime;

SC_W_trans_nc_mcps_indiv = SC_w_star .* dy_v .* dx_u;
SC_W_trans_nc_mcps = nansum(SC_W_trans_nc_mcps_indiv,1);
SC_W_trans_nc = SC_W_trans_nc_mcps * 10^-6;

SC_W_trans = NaN(1, length(lon_v));
SC_W_trans(57) = SC_W_trans_nc(57);
for jj = ALLC_lon_u_ind(2:end)
    SC_W_trans(jj) = ...
        SC_W_trans(jj-1) + SC_W_trans_nc(jj);
    
end


%
SC_W_trans_prime = NaN(1, length(lon_v));
for jj = ALLC_lon_u_ind(1:end-1)
    SC_W_trans_prime(jj) = - (-SC_U_trans(jj+1) + SC_U_trans(57) ...
        + SC_north_trans(jj) + SC_south_trans(jj));
    
end


% 
SC_U_trans_prime = NaN(1, length(lon_u));
SC_U_trans_prime(57) = SC_U_trans(57);
for jj = ALLC_lon_u_ind(1:end-1)
    SC_U_trans_prime(jj+1) = SC_U_trans_prime(57) + ...
        + SC_north_trans(jj) + SC_south_trans(jj) + ...
        SC_W_trans(jj);
    
end


%%
% close all
fig1 = figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1])

plot(lon_u, SC_U_trans, 'linewidth', 2)

hold all
plot(lon_u, SC_U_trans_prime, 'linewidth', 2)
plot(lon_v, SC_north_trans, '--', 'linewidth', 2)
plot(lon_v, SC_south_trans, '--', 'linewidth', 2)
plot(lon_v, SC_W_trans, ':', 'linewidth', 2)
plot(lon_v, SC_W_trans_prime, ':', 'linewidth', 2)


legend('U_{SBC}', 'U''_{SBC}', 'V_{SA}', 'V_{SBC}', 'W_{SBC}', 'W''_{SBC}')
axis([115 147 -2 2])
grid

% export_fig(fig1, ['Dropbox/SACS_work/figures/' ...
% 'm18_calc_transport/trans'], ...
%     '-m2', '-nocrop')
% 
% close

