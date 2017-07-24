%% get SBC and DRC
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/KDau_currents'])
load([data_path 'SACS_data/aus8_currents'])
load([data_path 'SACS_data/KDau_U_prime'])
load([data_path 'SACS_data/KDau_V_prime'])

a = aus8_coor.a;
pi180 = aus8_coor.pi180;
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;
U_prime_up = KDau_currents.ptop_to_pmid.U_prime.mean;
V_prime_up = KDau_currents.ptop_to_pmid.V_prime.mean;
U_g_prime_dw = KDau_currents.pmid_to_pbot.U_g_prime.mean;
V_g_prime_dw = KDau_currents.pmid_to_pbot.V_g_prime.mean;
lat_v_SBC_north = aus8_currents.lat_v_SBC_north;
lat_v_SBC_south = aus8_currents.lat_v_SBC_south;
lat_v_DRC_north = aus8_currents.lat_v_DRC_north;
lat_v_DRC_south = aus8_currents.lat_v_DRC_south;
lon_u_ALLC = aus8_currents.lon_u_ALLC;


%% Prep
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
KDau_U_prime.mean(isnan(KDau_U_prime.mean)) = 0;
KDau_V_prime.mean(isnan(KDau_V_prime.mean)) = 0;
du = KDau_U_prime.mean(:,2:end) - KDau_U_prime.mean(:,1:end-1);
dudx = du ./ dx_u;
lat_v_repmat = repmat(lat_v,1,length(lon_v));
dv = KDau_V_prime.mean(1:end-1,:).* ...
    cos(lat_v_repmat(1:end-1,:) * pi180) - ...
    KDau_V_prime.mean(2:end,:).* ...
    cos(lat_v_repmat(2:end,:) * pi180);
lat_u_repmat = repmat(lat_u,1,length(lon_v));
dvdy = 1./cos(lat_u_repmat * pi180).*dv./dy_v;
div_UV_prime.mean = dudx + dvdy;


%% SBC Ut
Ut_prime_up_SBC = zeros(size(Ut_prime_up));

% First
lat_north_ind_1 = ...
    find(ismember(lat_v, lat_v_SBC_north_repelem(1)));
lat_south_ind_1 = ...
    find(ismember(lat_v, lat_v_SBC_south_repelem(1)));
lat_vec_now = lat_north_ind_1:lat_south_ind_1-1;
lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(1));
Ut_prime_up_SBC(lat_vec_now, ...
    lon_u_ALLC_repelem_ind) = ...
    Ut_prime_up(lat_vec_now,...
    lon_u_ALLC_repelem_ind);
% Last
lat_north_ind_1 = ...
    find(ismember(lat_v, lat_v_SBC_north_repelem(end)));
lat_south_ind_1 = ...
    find(ismember(lat_v, lat_v_SBC_south_repelem(end)));
lat_vec_now = lat_north_ind_1:lat_south_ind_1-1;
lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(end));
Ut_prime_up_SBC(lat_vec_now, ...
    lon_u_ALLC_repelem_ind) = ...
    Ut_prime_up(lat_vec_now,...
    lon_u_ALLC_repelem_ind);

jj_count = 0;
for jj = 2 : 2 : length(lon_u_ALLC_repelem(1:end-2))
    jj_count = jj_count + 1;
    
    lat_north_ind_1 = ...
        find(ismember(lat_v, lat_v_SBC_north_repelem(jj)));
    lat_north_ind_2 = ...
        find(ismember(lat_v, lat_v_SBC_north_repelem(jj+1)));
    lat_north_12 = [lat_north_ind_1, lat_north_ind_2];
    lat_south_ind_1 = ...
        find(ismember(lat_v, lat_v_SBC_south_repelem(jj)));
    lat_south_ind_2 = ...
        find(ismember(lat_v, lat_v_SBC_south_repelem(jj+1)));
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
    lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(jj));
    Ut_prime_up_SBC(lat_vec_now, ...
        lon_u_ALLC_repelem_ind) = ...
        Ut_prime_up(lat_vec_now,...
        lon_u_ALLC_repelem_ind);
end

SBC_Ut_mcps = sum(Ut_prime_up_SBC, 1);
SBC_Ut = SBC_Ut_mcps .* 10^-6;


%% DRC Ut up
Ut_prime_up_DRC = zeros(size(Ut_prime_up));

% First
lat_north_ind_1 = ...
    find(ismember(lat_v, lat_v_DRC_north_repelem(1)));
lat_south_ind_1 = ...
    find(ismember(lat_v, lat_v_DRC_south_repelem(1)));
lat_vec_now = lat_north_ind_1:lat_south_ind_1-1;
lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(1));
Ut_prime_up_DRC(lat_vec_now, ...
    lon_u_ALLC_repelem_ind) = ...
    Ut_prime_up(lat_vec_now,...
    lon_u_ALLC_repelem_ind);
% Last
lat_north_ind_1 = ...
    find(ismember(lat_v, lat_v_DRC_north_repelem(end)));
lat_south_ind_1 = ...
    find(ismember(lat_v, lat_v_DRC_south_repelem(end)));
lat_vec_now = lat_north_ind_1:lat_south_ind_1-1;
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
        find(ismember(lat_v, lat_v_DRC_south_repelem(jj)));
    lat_south_ind_2 = ...
        find(ismember(lat_v, lat_v_DRC_south_repelem(jj-1)));
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
    lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(jj));
    Ut_prime_up_DRC(lat_vec_now, ...
        lon_u_ALLC_repelem_ind) = ...
        Ut_prime_up(lat_vec_now,...
        lon_u_ALLC_repelem_ind);
end

DRC_Ut_mcps_up = sum(Ut_prime_up_DRC, 1);
DRC_Ut_up = DRC_Ut_mcps_up .* 10^-6;


%% DRC Ut dw
Ut_g_prime_dw_DRC = zeros(size(Ut_g_prime_dw));

% First
lat_north_ind_1 = ...
    find(ismember(lat_v, lat_v_SBC_north_repelem(1)));
lat_south_ind_1 = ...
    find(ismember(lat_v, lat_v_DRC_south_repelem(1)));
lat_vec_now = lat_north_ind_1:lat_south_ind_1-1;
lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(1));
Ut_g_prime_dw_DRC(lat_vec_now, ...
    lon_u_ALLC_repelem_ind) = ...
    Ut_g_prime_dw(lat_vec_now,...
    lon_u_ALLC_repelem_ind);
% Last
lat_north_ind_1 = ...
    find(ismember(lat_v, lat_v_SBC_north_repelem(end)));
lat_south_ind_1 = ...
    find(ismember(lat_v, lat_v_DRC_south_repelem(end)));
lat_vec_now = lat_north_ind_1:lat_south_ind_1-1;
lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(end));
Ut_g_prime_dw_DRC(lat_vec_now, ...
    lon_u_ALLC_repelem_ind) = ...
    Ut_g_prime_dw(lat_vec_now,...
    lon_u_ALLC_repelem_ind);

jj_count = 0;
for jj = length(lon_u_ALLC_repelem(1:end-1)) : -2 : 3
    jj_count = jj_count + 1;
    
    lat_north_ind_1 = ...
        find(ismember(lat_v, lat_v_SBC_north_repelem(jj)));
    lat_north_ind_2 = ...
        find(ismember(lat_v, lat_v_SBC_north_repelem(jj-1)));
    lat_north_12 = [lat_north_ind_1, lat_north_ind_2];
    lat_south_ind_1 = ...
        find(ismember(lat_v, lat_v_DRC_south_repelem(jj)));
    lat_south_ind_2 = ...
        find(ismember(lat_v, lat_v_DRC_south_repelem(jj-1)));
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
    lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(jj));
    Ut_g_prime_dw_DRC(lat_vec_now, ...
        lon_u_ALLC_repelem_ind) = ...
        Ut_g_prime_dw(lat_vec_now,...
        lon_u_ALLC_repelem_ind);
end

DRC_Ut_mcps_dw = sum(Ut_g_prime_dw_DRC, 1);
DRC_Ut_dw = DRC_Ut_mcps_dw .* 10^-6;


%% DRC Ut
DRC_Ut = DRC_Ut_up + DRC_Ut_dw;


%% SBC Vtn
Vtn_prime_up_SBC = zeros(1, length(lon_v));

% do the first
lon_first_idx = 57;
lat_ind_V_Vtn_prime_up_SBC = ...
    lat_v == lat_v_SBC_north(1);
V_Vtn_prime_up_SBC_now = ...
    Vt_prime_up(lat_ind_V_Vtn_prime_up_SBC, lon_first_idx);
V_Vtn_prime_up_SBC = -V_Vtn_prime_up_SBC_now;
Vtn_prime_up_SBC(lon_first_idx) = V_Vtn_prime_up_SBC;

% and last cases
lon_last_idx = 312;
lat_ind_V_Vtn_prime_up_SBC = ...
    lat_v == lat_v_SBC_north(end);
V_Vtn_prime_up_SBC_now = ...
    Vt_prime_up(lat_ind_V_Vtn_prime_up_SBC, lon_last_idx);
V_Vtn_prime_up_SBC = -V_Vtn_prime_up_SBC_now;
Vtn_prime_up_SBC(lon_last_idx) = V_Vtn_prime_up_SBC;

% in between cases
for jj = 2 : length(lon_u_ALLC)-2
    lon_u_ind = find(lon_u == lon_u_ALLC(jj));
    
    %
    lat_ind_V_Vtn_prime_up_SBC = lat_v == lat_v_SBC_north(jj);
    V_Vtn_prime_up_SBC_now = ...
        Vt_prime_up(lat_ind_V_Vtn_prime_up_SBC, lon_u_ind);
    V_Vtn_prime_up_SBC = -V_Vtn_prime_up_SBC_now;
    
    %
    [U_Vtn_prime_up_SBC_prev, U_Vtn_prime_up_SBC_next] = deal(0);
    
    % if the previous lat is to the south of the current one
    if lat_v_SBC_north(jj-1) < lat_v_SBC_north(jj)
        %
        lat_u_within_ind = find(...
            lat_u < lat_v_SBC_north(jj) & ...
            lat_u > lat_v_SBC_north(jj-1));
        U_Vtn_prime_up_SBC_all_now = ...
            Ut_prime_up(lat_u_within_ind, lon_u_ind);
        U_Vtn_prime_up_SBC_all = U_Vtn_prime_up_SBC_all_now;
        U_Vtn_prime_up_SBC_prev = nansum(U_Vtn_prime_up_SBC_all);
        if isnan(U_Vtn_prime_up_SBC_prev)
            U_Vtn_prime_up_SBC_prev = 0;
        end
    end
    
    % if the next lat is to the south of the current one
    if lat_v_SBC_north(jj+1) < lat_v_SBC_north(jj)
        %
        lat_u_within_ind = find(...
            lat_u < lat_v_SBC_north(jj) & ...
            lat_u > lat_v_SBC_north(jj+1));
        U_Vtn_prime_up_SBC_all_now = ...
            Ut_prime_up(lat_u_within_ind, lon_u_ind+1);
        U_Vtn_prime_up_SBC_all = -U_Vtn_prime_up_SBC_all_now;
        U_Vtn_prime_up_SBC_next = nansum(U_Vtn_prime_up_SBC_all);
    end
    if isnan(U_Vtn_prime_up_SBC_next)
        U_Vtn_prime_up_SBC_next = 0;
    end
    
    %
    Vtn_prime_up_SBC(lon_u_ind) = V_Vtn_prime_up_SBC + ...
        U_Vtn_prime_up_SBC_prev + U_Vtn_prime_up_SBC_next;
end


%% SBC Vts
Vts_prime_up_SBC = zeros(1, length(lon_v));

% do the first
lon_first_idx = 57;
lat_ind_V_Vts_prime_up_SBC = ...
    lat_v == lat_v_SBC_south(1);
V_Vts_prime_up_SBC_now = ...
    Vt_prime_up(lat_ind_V_Vts_prime_up_SBC, lon_first_idx);
V_Vts_prime_up_SBC = V_Vts_prime_up_SBC_now;
Vts_prime_up_SBC(lon_first_idx) = V_Vts_prime_up_SBC;

% and last cases
lon_last_idx = 312;
lat_ind_V_Vts_prime_up_SBC = ...
    lat_v == lat_v_SBC_south(end);
V_Vts_prime_up_SBC_now = ...
    Vt_prime_up(lat_ind_V_Vts_prime_up_SBC, lon_last_idx);
V_Vts_prime_up_SBC = V_Vts_prime_up_SBC_now;
Vts_prime_up_SBC(lon_last_idx) = V_Vts_prime_up_SBC;

% in between cases
for jj = 2 : length(lon_u_ALLC)-2
    lon_u_ind = find(lon_u == lon_u_ALLC(jj));
    
    %
    lat_ind_V_Vts_prime_up_SBC = lat_v == lat_v_SBC_south(jj);
    V_Vts_prime_up_SBC_now = ...
        Vt_prime_up(lat_ind_V_Vts_prime_up_SBC, lon_u_ind);
    V_Vts_prime_up_SBC = V_Vts_prime_up_SBC_now;
    
    %
    [U_Vts_prime_up_SBC_prev, U_Vts_prime_up_SBC_next] = deal(0);
    
    % if the previous lat is to the north of the current one
    if lat_v_SBC_south(jj-1) > lat_v_SBC_south(jj)
        %
        lat_u_within_ind = find(...
            lat_u > lat_v_SBC_south(jj) & ...
            lat_u < lat_v_SBC_south(jj-1));
        U_Vts_prime_up_SBC_all_now = ...
            Ut_prime_up(lat_u_within_ind, lon_u_ind);
        U_Vts_prime_up_SBC_all = U_Vts_prime_up_SBC_all_now;
        U_Vts_prime_up_SBC_prev = nansum(U_Vts_prime_up_SBC_all);
    end
    
    % if the next lat is to the north of the current one
    if lat_v_SBC_south(jj+1) > lat_v_SBC_south(jj)
        %
        lat_u_within_ind = find(...
            lat_u > lat_v_SBC_south(jj) & ...
            lat_u < lat_v_SBC_south(jj+1));
        U_Vts_prime_up_SBC_all_now = ...
            Ut_prime_up(lat_u_within_ind, lon_u_ind+1);
        U_Vts_prime_up_SBC_all = -U_Vts_prime_up_SBC_all_now;
        U_Vts_prime_up_SBC_next = nansum(U_Vts_prime_up_SBC_all);
    end
    
    %
    Vts_prime_up_SBC(lon_u_ind) = ...
        V_Vts_prime_up_SBC + U_Vts_prime_up_SBC_prev + ...
        U_Vts_prime_up_SBC_next;
end


%% SBC Vtc
SBC_Vtn = Vtn_prime_up_SBC .* 10^-6;
SBC_Vts = Vts_prime_up_SBC .* 10^-6;

[SBC_Vtnc, SBC_Vtsc] = deal(zeros(1, length(lon_v)));
SBC_Vtnc(57) = SBC_Vtn(57);
SBC_Vtsc(57) = SBC_Vts(57);

for jj = 58 : 312
    SBC_Vtnc(jj) = SBC_Vtnc(jj-1) + SBC_Vtn(jj);
    SBC_Vtsc(jj) = SBC_Vtsc(jj-1) + SBC_Vts(jj);
end


%% DRC Vtn up
Vtn_prime_up_DRC = zeros(1, length(lon_v));

% do the first
lon_first_idx = 57;
lat_ind_V_Vtn_prime_up_DRC = ...
    lat_v == lat_v_DRC_north(1);
V_Vtn_prime_up_DRC_now = ...
    Vt_prime_up(lat_ind_V_Vtn_prime_up_DRC, lon_first_idx);
V_Vtn_prime_up_DRC = -V_Vtn_prime_up_DRC_now;
Vtn_prime_up_DRC(lon_first_idx) = V_Vtn_prime_up_DRC;

% and last cases
lon_last_idx = 312;
lat_ind_V_Vtn_prime_up_DRC = ...
    lat_v == lat_v_DRC_north(end);
V_Vtn_prime_up_DRC_now = ...
    Vt_prime_up(lat_ind_V_Vtn_prime_up_DRC, lon_last_idx);
V_Vtn_prime_up_DRC = -V_Vtn_prime_up_DRC_now;
Vtn_prime_up_DRC(lon_last_idx) = V_Vtn_prime_up_DRC;

% in between cases
for jj = length(lon_u_ALLC)-2 : -1 : 2
    lon_u_ind = find(lon_u == lon_u_ALLC(jj));
    
    %
    lat_ind_V_Vtn_prime_up_DRC = lat_v == lat_v_DRC_north(jj);
    V_Vtn_prime_up_DRC_now = ...
        Vt_prime_up(lat_ind_V_Vtn_prime_up_DRC, lon_u_ind);
    V_Vtn_prime_up_DRC = -V_Vtn_prime_up_DRC_now;
    
    %
    [U_Vtn_prime_up_DRC_prev, U_Vtn_prime_up_DRC_next] = deal(0);
    
    % if the previous lat is to the south of the current one
    if lat_v_DRC_north(jj+1) < lat_v_DRC_north(jj)
        %
        lat_u_within_ind = find(...
            lat_u < lat_v_DRC_north(jj) & ...
            lat_u > lat_v_DRC_north(jj+1));
        U_Vtn_prime_up_DRC_all_now = ...
            Ut_prime_up(lat_u_within_ind, lon_u_ind);
        U_Vtn_prime_up_DRC_all = -U_Vtn_prime_up_DRC_all_now;
        U_Vtn_prime_up_DRC_prev = nansum(U_Vtn_prime_up_DRC_all);
        if isnan(U_Vtn_prime_up_DRC_prev)
            U_Vtn_prime_up_DRC_prev = 0;
        end
    end
    
    % if the next lat is to the south of the current one
    if lat_v_SBC_north(jj-1) < lat_v_SBC_north(jj)
        %
        lat_u_within_ind = find(...
            lat_u < lat_v_DRC_north(jj) & ...
            lat_u > lat_v_DRC_north(jj-1));
        U_Vtn_prime_up_DRC_all_now = ...
            Ut_prime_up(lat_u_within_ind, lon_u_ind+1);
        U_Vtn_prime_up_DRC_all = U_Vtn_prime_up_DRC_all_now;
        U_Vtn_prime_up_DRC_next = nansum(U_Vtn_prime_up_DRC_all);
    end
    if isnan(U_Vtn_prime_up_DRC_next)
        U_Vtn_prime_up_DRC_next = 0;
    end
    
    %
    Vtn_prime_up_DRC(lon_u_ind) = V_Vtn_prime_up_DRC + ...
        U_Vtn_prime_up_DRC_prev + U_Vtn_prime_up_DRC_next;
end


%% DRC Vtn dw
% % This should be zero...
% Vtn_g_prime_dw_DRC = zeros(1, length(lon_v));

% % do the first
% lon_first_idx = 57;
% lat_ind_V_Vtn_g_prime_dw_DRC = ...
%     lat_v == lat_v_SBC_north(1);
% V_Vtn_g_prime_dw_DRC_now = ...
%     Vt_g_prime_dw(lat_ind_V_Vtn_g_prime_dw_DRC, lon_first_idx);
% V_Vtn_g_prime_dw_DRC = -V_Vtn_g_prime_dw_DRC_now;
% Vtn_g_prime_dw_DRC(lon_first_idx) = V_Vtn_g_prime_dw_DRC;
% 
% % and last cases
% lon_last_idx = 312;
% lat_ind_V_Vtn_g_prime_dw_DRC = ...
%     lat_v == lat_v_SBC_north(end);
% V_Vtn_g_prime_dw_DRC_now = ...
%     Vt_g_prime_dw(lat_ind_V_Vtn_g_prime_dw_DRC, lon_last_idx);
% V_Vtn_g_prime_dw_DRC = -V_Vtn_g_prime_dw_DRC_now;
% Vtn_g_prime_dw_DRC(lon_last_idx) = V_Vtn_g_prime_dw_DRC;
% 
% % in between cases
% for jj = 2 : length(lon_u_ALLC)-2
%     lon_u_ind = find(lon_u == lon_u_ALLC(jj));
%     
%     %
%     lat_ind_V_Vtn_g_prime_dw_DRC = lat_v == lat_v_SBC_north(jj);
%     V_Vtn_g_prime_dw_DRC_now = ...
%         Vt_g_prime_dw(lat_ind_V_Vtn_g_prime_dw_DRC, lon_u_ind);
%     V_Vtn_g_prime_dw_DRC = -V_Vtn_g_prime_dw_DRC_now;
%     
%     %
%     [U_Vtn_g_prime_dw_DRC_prev, U_Vtn_g_prime_dw_DRC_next] = deal(0);
%     
%     % if the previous lat is to the south of the current one
%     if lat_v_SBC_north(jj-1) < lat_v_SBC_north(jj)
%         %
%         lat_u_within_ind = find(...
%             lat_u < lat_v_SBC_north(jj) & ...
%             lat_u > lat_v_SBC_north(jj-1));
%         U_Vtn_g_prime_dw_DRC_all_now = ...
%             Ut_g_prime_dw(lat_u_within_ind, lon_u_ind);
%         U_Vtn_g_prime_dw_DRC_all = U_Vtn_g_prime_dw_DRC_all_now;
%         U_Vtn_g_prime_dw_DRC_prev = nansum(U_Vtn_g_prime_dw_DRC_all);
%         if isnan(U_Vtn_g_prime_dw_DRC_prev)
%             U_Vtn_g_prime_dw_DRC_prev = 0;
%         end
%     end
%     
%     % if the next lat is to the south of the current one
%     if lat_v_SBC_north(jj+1) < lat_v_SBC_north(jj)
%         %
%         lat_u_within_ind = find(...
%             lat_u < lat_v_SBC_north(jj) & ...
%             lat_u > lat_v_SBC_north(jj+1));
%         U_Vtn_g_prime_dw_DRC_all_now = ...
%             Ut_g_prime_dw(lat_u_within_ind, lon_u_ind+1);
%         U_Vtn_g_prime_dw_DRC_all = -U_Vtn_g_prime_dw_DRC_all_now;
%         U_Vtn_g_prime_dw_DRC_next = nansum(U_Vtn_g_prime_dw_DRC_all);
%     end
%     if isnan(U_Vtn_g_prime_dw_DRC_next)
%         U_Vtn_g_prime_dw_DRC_next = 0;
%     end
%     
%     %
%     Vtn_g_prime_dw_DRC(lon_u_ind) = V_Vtn_g_prime_dw_DRC + ...
%         U_Vtn_g_prime_dw_DRC_prev + U_Vtn_g_prime_dw_DRC_next;
% end
% 
% for n = 1 : length(lon_u)-1
%     if Vtn_g_prime_dw_DRC(n) ~= 0
%         error('Vtn_g_prime_dw_DRC = 0 somewhere :( !')
%     end
% end
% disp('Vtn_g_prime_dw_DRC = 0 everywhere :)')


%% DRC Vts up
Vts_prime_up_DRC = zeros(1, length(lon_v));

% do the first
lon_first_idx = 57;
lat_ind_V_Vts_prime_up_DRC = ...
    lat_v == lat_v_DRC_south(1);
V_Vts_prime_up_DRC_now = ...
    Vt_prime_up(lat_ind_V_Vts_prime_up_DRC, lon_first_idx);
V_Vts_prime_up_DRC = V_Vts_prime_up_DRC_now;
Vts_prime_up_DRC(lon_first_idx) = V_Vts_prime_up_DRC;

% and last cases
lon_last_idx = 312;
lat_ind_V_Vts_prime_up_DRC = ...
    lat_v == lat_v_DRC_south(end);
V_Vts_prime_up_DRC_now = ...
    Vt_prime_up(lat_ind_V_Vts_prime_up_DRC, lon_last_idx);
V_Vts_prime_up_DRC = V_Vts_prime_up_DRC_now;
Vts_prime_up_DRC(lon_last_idx) = V_Vts_prime_up_DRC;

% in between cases
for jj = length(lon_u_ALLC)-2 : -1 : 2
    lon_u_ind = find(lon_u == lon_u_ALLC(jj));
    
    %
    lat_ind_V_Vts_prime_up_DRC = lat_v == lat_v_DRC_south(jj);
    V_Vts_prime_up_DRC_now = ...
        Vt_prime_up(lat_ind_V_Vts_prime_up_DRC, lon_u_ind);
    V_Vts_prime_up_DRC = V_Vts_prime_up_DRC_now;
    
    %
    [U_Vts_prime_up_DRC_prev, U_Vts_prime_up_DRC_next] = deal(0);
    
    % if the previous lat is to the north of the current one
    if lat_v_DRC_south(jj+1) > lat_v_DRC_south(jj)
        %
        lat_u_within_ind = find(...
            lat_u > lat_v_DRC_south(jj) & ...
            lat_u < lat_v_DRC_south(jj+1));
        U_Vts_prime_up_DRC_all_now = ...
            Ut_prime_up(lat_u_within_ind, lon_u_ind+1);
        U_Vts_prime_up_DRC_all = -U_Vts_prime_up_DRC_all_now;
        U_Vts_prime_up_DRC_prev = nansum(U_Vts_prime_up_DRC_all);
    end
    
    % if the next lat is to the north of the current one
    if lat_v_DRC_south(jj-1) > lat_v_DRC_south(jj)
        %
        lat_u_within_ind = find(...
            lat_u > lat_v_DRC_south(jj) & ...
            lat_u < lat_v_DRC_south(jj-1));
        U_Vts_prime_up_DRC_all_now = ...
            Ut_prime_up(lat_u_within_ind, lon_u_ind);
        U_Vts_prime_up_DRC_all = U_Vts_prime_up_DRC_all_now;
        U_Vts_prime_up_DRC_next = nansum(U_Vts_prime_up_DRC_all);
    end
    
    %
    Vts_prime_up_DRC(lon_u_ind) = ...
        V_Vts_prime_up_DRC + ...
        U_Vts_prime_up_DRC_prev + U_Vts_prime_up_DRC_next;
end


%% DRC Vts dw
Vts_g_prime_dw_DRC = zeros(1, length(lon_v));

% do the first
lon_first_idx = 57;
lat_ind_V_Vts_g_prime_dw_DRC = ...
    lat_v == lat_v_DRC_south(1);
V_Vts_g_prime_dw_DRC_now = ...
    Vt_g_prime_dw(lat_ind_V_Vts_g_prime_dw_DRC, lon_first_idx);
V_Vts_g_prime_dw_DRC = V_Vts_g_prime_dw_DRC_now;
Vts_g_prime_dw_DRC(lon_first_idx) = V_Vts_g_prime_dw_DRC;

% and last cases
lon_last_idx = 312;
lat_ind_V_Vts_g_prime_dw_DRC = ...
    lat_v == lat_v_DRC_south(end);
V_Vts_g_prime_dw_DRC_now = ...
    Vt_g_prime_dw(lat_ind_V_Vts_g_prime_dw_DRC, lon_last_idx);
V_Vts_g_prime_dw_DRC = V_Vts_g_prime_dw_DRC_now;
Vts_g_prime_dw_DRC(lon_last_idx) = V_Vts_g_prime_dw_DRC;

% in between cases
for jj = length(lon_u_ALLC)-2 : -1 : 2
    lon_u_ind = find(lon_u == lon_u_ALLC(jj));
    
    %
    lat_ind_V_Vts_g_prime_dw_DRC = lat_v == lat_v_DRC_south(jj);
    V_Vts_g_prime_dw_DRC_now = ...
        Vt_g_prime_dw(lat_ind_V_Vts_g_prime_dw_DRC, lon_u_ind);
    V_Vts_g_prime_dw_DRC = V_Vts_g_prime_dw_DRC_now;
    
    %
    [U_Vts_g_prime_dw_DRC_prev, U_Vts_g_prime_dw_DRC_next] = deal(0);
    
    % if the previous lat is to the north of the current one
    if lat_v_DRC_south(jj+1) > lat_v_DRC_south(jj)
        %
        lat_u_within_ind = find(...
            lat_u > lat_v_DRC_south(jj) & ...
            lat_u < lat_v_DRC_south(jj+1));
        U_Vts_g_prime_dw_DRC_all_now = ...
            Ut_g_prime_dw(lat_u_within_ind, lon_u_ind+1);
        U_Vts_g_prime_dw_DRC_all = -U_Vts_g_prime_dw_DRC_all_now;
        U_Vts_g_prime_dw_DRC_prev = nansum(U_Vts_g_prime_dw_DRC_all);
    end
    
    % if the next lat is to the north of the current one
    if lat_v_DRC_south(jj-1) > lat_v_DRC_south(jj)
        %
        lat_u_within_ind = find(...
            lat_u > lat_v_DRC_south(jj) & ...
            lat_u < lat_v_DRC_south(jj-1));
        U_Vts_g_prime_dw_DRC_all_now = ...
            Ut_g_prime_dw(lat_u_within_ind, lon_u_ind);
        U_Vts_g_prime_dw_DRC_all = U_Vts_g_prime_dw_DRC_all_now;
        U_Vts_g_prime_dw_DRC_next = nansum(U_Vts_g_prime_dw_DRC_all);
    end
    
    %
    Vts_g_prime_dw_DRC(lon_u_ind) = ...
        V_Vts_g_prime_dw_DRC + ...
        U_Vts_g_prime_dw_DRC_prev + U_Vts_g_prime_dw_DRC_next;
end


%% DRC Vtc
DRC_Vtn = Vtn_prime_up_DRC .* 10^-6;
DRC_Vts_up = Vts_prime_up_DRC .* 10^-6;
DRC_Vts_dw = Vts_g_prime_dw_DRC .* 10^-6;

[DRC_Vtnc, DRC_Vtsc_up, DRC_Vtsc_dw] = deal(zeros(1, length(lon_v)));
DRC_Vtnc(312) = DRC_Vtn(312);
DRC_Vtsc_up(312) = DRC_Vts_up(312);
DRC_Vtsc_dw(312) = DRC_Vts_dw(312);

for jj = 311 : -1 : 57
    DRC_Vtnc(jj) = DRC_Vtnc(jj+1) + DRC_Vtn(jj);
    DRC_Vtsc_up(jj) = DRC_Vtsc_up(jj+1) + DRC_Vts_up(jj);
    DRC_Vtsc_dw(jj) = DRC_Vtsc_dw(jj+1) + DRC_Vts_dw(jj);
end


%% DRC Vts
DRC_Vts = DRC_Vts_up + DRC_Vts_dw;
DRC_Vtsc = DRC_Vtsc_up + DRC_Vtsc_dw;


%% SBC Wt
% Get Ut_prime_up_SBC Vt_prime_up_SBC for W
V_prime_up_SBC_for_W = NaN(size(V_prime_up));
U_prime_up_SBC_for_W = NaN(size(U_prime_up));

nc = 0;
for n = 57 : 312
    nc = nc + 1;
    lat_v_SBC_north_ind = find(lat_v==lat_v_SBC_north(nc));
    lat_v_SBC_south_ind = find(lat_v==lat_v_SBC_south(nc));
    lat_v_SBC_ind = lat_v_SBC_north_ind : lat_v_SBC_south_ind;
    
    lat_u_SBC_north_ind = find(lat_u<lat_v_SBC_north(nc), 1, 'first');
    lat_u_SBC_south_ind = find(lat_u>lat_v_SBC_south(nc), 1, 'last');
    lat_u_SBC_ind = lat_u_SBC_north_ind : lat_u_SBC_south_ind;
    
    V_prime_up_SBC_for_W(lat_v_SBC_ind,n) = ...
        V_prime_up(lat_v_SBC_ind,n);
    U_prime_up_SBC_for_W(lat_u_SBC_ind,n:n+1) = ...
        U_prime_up(lat_u_SBC_ind,n:n+1);
end

du = U_prime_up_SBC_for_W(:,2:end) - ...
    U_prime_up_SBC_for_W(:,1:end-1);
dudx = du ./ dx_u;
dv = V_prime_up_SBC_for_W(1:end-1,:).* ...
    cos(lat_v_repmat(1:end-1,:) * pi180) - ...
    V_prime_up_SBC_for_W(2:end,:).* ...
    cos(lat_v_repmat(2:end,:) * pi180);
dvdy = 1./cos(lat_u_repmat * pi180).*dv./dy_v;
div_UV_prime_up_SBC = dudx + dvdy;

% for ii = 1 : length(lat_u)
%     for jj = 1 : length(lon_v)
%         if isfinite(div_UV_prime_up_SBC(ii,jj))
%             if div_UV_prime_up_SBC(ii,jj) == div_UV_prime.mean(ii,jj)
%                 div_UV_prime_up_SBC(ii,jj) = 0;
%             end
%         end
%     end
% end
% disp('div_UV_prime_up_SBC was corrected')

w_prime_up_SBC = div_UV_prime_up_SBC;
W_prime_up_mcps_SBC = w_prime_up_SBC .* dy_v .* dx_u;
SBC_Wt_mcps = nansum(W_prime_up_mcps_SBC,1);
SBC_Wt = SBC_Wt_mcps * 10^-6;


%% DRC Wt
DRC_Wt = -SBC_Wt;


%% Wtc
[SBC_Wtc, DRC_Wtc] = deal(zeros(1, length(lon_v)));

SBC_Wtc(57) = SBC_Wt(57);
for jj = 58 : 312
    SBC_Wtc(jj) = SBC_Wtc(jj-1) + SBC_Wt(jj);
end
DRC_Wtc(312) = DRC_Wt(312);
for jj = 311 : -1 : 57
    DRC_Wtc(jj) = DRC_Wtc(jj+1) + DRC_Wt(jj);
end


%% SBC Ut star
SBC_Ut_star = zeros(1, length(lon_u));
SBC_Ut_star(57) = SBC_Ut(57);
for jj = 57 : 312
    SBC_Ut_star(jj+1) = SBC_Ut_star(57) + ...
        + SBC_Vtnc(jj) + SBC_Vtsc(jj) + SBC_Wtc(jj);
end


%% DRC Ut star
DRC_Ut_star = zeros(1, length(lon_u));
DRC_Ut_star(313) = DRC_Ut(313);
for jj = 312 : -1 : 57
    DRC_Ut_star(jj) = DRC_Ut_star(313) + ...
        -(DRC_Vtnc(jj) + DRC_Vtsc(jj) + DRC_Wtc(jj));
end


%%
close all
fig_n = 1;
rowcols = [1 2];
rowcols_size = [20 15]; % cm
margs = [1 1 1 1]; % cm
gaps = [1.5 1]; % cm
cmaps_levels = 12;

font_size = 12;
fig_color = [1 1 1];

fig = figure(fig_n);
rowN = rowcols(1); colN = rowcols(2);
[rm, cm] = meshgrid(rowN:-1:1, 1:colN);
x_sp = rowcols_size(1); % x subplot length
y_sp = rowcols_size(2); % y subplot length
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = gaps(1); % gap width between subplots
gap_h = gaps(2); % gap height between subplots
marg_b = margs(3); % bottom_KDau margin
marg_t = margs(4); % top margin
marg_l = margs(1); % left margin
marg_r = margs(2); % right margin
set(fig,'units','centimeters',...
    'position',[0 0 ...
    (marg_l+colN*x_sp+gap_w*(colN-1)+marg_r) ...
    (marg_b+rowN*y_sp+gap_h*(rowN-1)+marg_t)]);

sp = 1;
axes('Units','centimeters', ...
    'Position',[...
    (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)), ...
    (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
    x_sp, ...
    y_sp])
plot(lon_u, SBC_Ut, 'b-', 'linewidth', 1)
hold all
plot(lon_u, SBC_Ut_star, 'b--', 'linewidth', 1)
plot(lon_v, SBC_Vtnc, 'g-', 'linewidth', 1)
plot(lon_v, SBC_Vtsc, 'r-', 'linewidth', 1)
plot(lon_v, SBC_Wtc, 'm-', 'linewidth', 1)
legend(...
    'SBC $U_{t}$', ...
    'SBC $U_{t}^{*}$', ...
    'SBC $V_{t}nc$', ...
    'SBC $V_{t}sc$', ...
    'SBC $W_{t}c$')
axis([115 147 -3 3])
title('KDS75 SBC transports ($Sv$) VS longitude')
grid
set(gca,'layer','top','color',fig_color,...
    'fontsize',font_size,'tickdir','out')

sp = 2;
axes('Units','centimeters', ...
    'Position',[...
    (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)), ...
    (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
    x_sp, ...
    y_sp])
plot(lon_u, DRC_Ut, 'b-', 'linewidth', 1)
hold all
plot(lon_u, DRC_Ut_star, 'b--', 'linewidth', 1)
plot(lon_v, DRC_Vtnc, 'r-', 'linewidth', 1)
plot(lon_v, DRC_Vtsc, 'k-', 'linewidth', 1)
plot(lon_v, DRC_Vtsc_up, 'k:', 'linewidth', 1)
plot(lon_v, DRC_Vtsc_dw, 'k--', 'linewidth', 1)
plot(lon_v, DRC_Wtc, 'm-', 'linewidth', 1)
legend(...
    'DRC $U_{t}$', ...
    'DRC $U_{t}^{*}$', ...
    'DRC $V_{t}nc$', ...
    'DRC $V_{t}sc$', ...
    'DRC $V_{t}sc$ up', ...
    'DRC $V_{t}sc$ dw', ...
    'DRC $W_{t}c$')
axis([115 147 -17 17])
title('KDS75 DRC transports ($Sv$) VS longitude')
grid
set(gca,'layer','top','color',fig_color,...
    'fontsize',font_size,'tickdir','out')

% outputls = ls(figures_path);
% scriptname = mfilename;
% if ~contains(outputls, scriptname)
%     mkdir(figures_path, scriptname)
% end
% export_fig(fig, ...
%     [figures_path mfilename '/' scriptname(1:3) ...
%     '_fig' num2str(fig_n) '_'], ...
%     '-m3')
% close


%%
close all
fig_n = 2;
rowcols = [1 2];
rowcols_size = [15 12.5]; % cm
margs = [1 1 1 1]; % cm
gaps = [1.5 1]; % cm
cmaps_levels = 12;

font_size = 12;
fig_color = [1 1 1];

fig = figure(fig_n);
rowN = rowcols(1); colN = rowcols(2);
[rm, cm] = meshgrid(rowN:-1:1, 1:colN);
x_sp = rowcols_size(1); % x subplot length
y_sp = rowcols_size(2); % y subplot length
col_ind = (repmat(1:colN,rowN,1))';
row_ind = (fliplr(repmat((1:rowN)',1,colN)))';
gap_w = gaps(1); % gap width between subplots
gap_h = gaps(2); % gap height between subplots
marg_b = margs(3); % bottom_KDau margin
marg_t = margs(4); % top margin
marg_l = margs(1); % left margin
marg_r = margs(2); % right margin
set(fig,'units','centimeters',...
    'position',[0 0 ...
    (marg_l+colN*x_sp+gap_w*(colN-1)+marg_r) ...
    (marg_b+rowN*y_sp+gap_h*(rowN-1)+marg_t)]);

sp = 1;
axes('Units','centimeters', ...
    'Position',[...
    (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)), ...
    (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
    x_sp, ...
    y_sp])
plot(lon_u, SBC_Ut-SBC_Ut_star, 'b-', 'linewidth', 1)
title('KDS75 SBC "leak": $U_{t}-U_{t}^{*}$ ($Sv$) VS longitude')
grid
set(gca,'layer','top','color',fig_color,...
    'fontsize',font_size,'tickdir','out')

sp = 2;
axes('Units','centimeters', ...
    'Position',[...
    (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)), ...
    (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
    x_sp, ...
    y_sp])
plot(lon_u, DRC_Ut-DRC_Ut_star, 'b-', 'linewidth', 1)
title('KDS75 DRC "leak": $U_{t}-U_{t}^{*}$ ($Sv$) VS longitude')
grid
set(gca,'layer','top','color',fig_color,...
    'fontsize',font_size,'tickdir','out')

% outputls = ls(figures_path);
% scriptname = mfilename;
% if ~contains(outputls, scriptname)
%     mkdir(figures_path, scriptname)
% end
% export_fig(fig, ...
%     [figures_path mfilename '/' scriptname(1:3) ...
%     '_fig' num2str(fig_n) '_'], ...
%     '-m3')
% close

