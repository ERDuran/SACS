%% get SBC and DRC
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_currents'])
load([data_path 'SACS_data/KDau_fulu'])
load([data_path 'SACS_data/KDau_fulv'])
load([data_path 'SACS_data/KDau_fulw'])

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
MTH = aus8_coor.MTH;
MTH{5} = 'mean';
depth = aus8_coor.depth;
depth_mid = aus8_coor.depth_mid;
depth_thkn = aus8_coor.depth_thkn;


%% 1) Define mid pressure and west, east start/end
%%% BOTTOM PRES
% pressure of the bottom surface currents
z_top = aus8_currents.z_top;
z_mid = aus8_currents.z_mid;
z_bot = aus8_currents.z_bot;
% index vector of the bottom surface currents
z_top_below_ind = find(depth==z_top);
z_mid_above_ind = find(depth==z_mid)-1;
z_mid_below_ind = find(depth==z_mid);
% pressure above the bottom
z_top_below = depth_mid(z_top_below_ind);
z_mid_above = depth_mid(z_mid_above_ind);
z_mid_below = depth_mid(z_mid_below_ind);
% index vector of currents pressure from surface to interface
z_mid_above_all_ind = depth_mid >= z_mid_above;
z_mid_below_all_ind = depth_mid <= z_mid_below;
%%% BOTTOM PRES


%% Prep
for t = 1 : 5    
    %% 3) Calculate UV_upper and UV_lower
    depth_thkn_perm = permute(depth_thkn, [3 2 1]);
    depth_thkn_u = repmat(depth_thkn_perm, [length(lat_u), length(lon_u)]);
    depth_thkn_v = repmat(depth_thkn_perm, [length(lat_v), length(lon_v)]);
    
    u_times_depth_thkn_u = KDau_fulu.(MTH{t}) .* depth_thkn_u;
    U_ztop_to_zmid = ...
        nansum(u_times_depth_thkn_u(:,:,z_mid_above_all_ind), 3);
    
    v_times_depth_thkn_v = KDau_fulv.(MTH{t}) .* depth_thkn_v;
    V_ztop_to_zmid = ...
        nansum(v_times_depth_thkn_v(:,:,z_mid_above_all_ind), 3);
    
    U_zmid_to_zbot = ...
        nansum(u_times_depth_thkn_u(:,:,z_mid_below_all_ind), 3);
    
    V_zmid_to_zbot = ...
        nansum(v_times_depth_thkn_v(:,:,z_mid_below_all_ind), 3);
    
    W_zmid = KDau_fulw.(MTH{t})(:,:,z_mid_above_ind+1);
    W_zbot = KDau_fulw.(MTH{t})(:,:,end);
    
    U_ztop_to_zmid(U_ztop_to_zmid==0) = NaN;
    V_ztop_to_zmid(V_ztop_to_zmid==0) = NaN;
    U_zmid_to_zbot(U_zmid_to_zbot==0) = NaN;
    V_zmid_to_zbot(V_zmid_to_zbot==0) = NaN;
    
    KDau_fcrt.ztop_to_zmid.fulu.(MTH{t}) = U_ztop_to_zmid;
    KDau_fcrt.ztop_to_zmid.fulv.(MTH{t}) = V_ztop_to_zmid;
    KDau_fcrt.zmid_to_zbot.fulu.(MTH{t}) = U_zmid_to_zbot;
    KDau_fcrt.zmid_to_zbot.fulv.(MTH{t}) = V_zmid_to_zbot;
    KDau_fcrt.zmid.fulw.(MTH{t}) = W_zmid;
    KDau_fcrt.zbot.fulw.(MTH{t}) = W_zbot;
    disp([MTH{t} ' OK!'])
    
    %
    fulU_up = U_ztop_to_zmid;
    fulV_up = V_ztop_to_zmid;
    fulU_dw = U_zmid_to_zbot;
    fulV_dw = V_zmid_to_zbot;
    fulW_mi = W_zmid;
    fulW_bo = W_zbot;
    
    %
    dx_v = NaN(length(lat_v), length(lon_u)-1);
    for ii = 1 : length(lat_v)
        dx_v(ii,:) = a * cos(lat_v(ii) * pi180) * ...
            (lon_u(2:end) - lon_u(1:end-1)) * pi180;
    end
    dy_raw = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
    dy_u = repmat(dy_raw, [1 length(lon_u)]);
    
    %
    dx_u = NaN(length(lat_u), length(lon_u)-1);
    for ii = 1 : length(lat_u)
        dx_u(ii,:) = a * cos(lat_u(ii) * pi180) * ...
            (lon_u(2:end) - lon_u(1:end-1)) * pi180;
    end
    dy_raw = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
    dy_v = repmat(dy_raw, [1 length(lon_v)]);
    
    %
    fulU_up(isnan(fulU_up)) = 0;
    fulV_up(isnan(fulV_up)) = 0;
    fulU_dw(isnan(fulU_dw)) = 0;
    fulV_dw(isnan(fulV_dw)) = 0;
    fulW_mi(isnan(fulW_mi)) = 0;
    fulW_bo(isnan(fulW_bo)) = 0;
    Ut_up = fulU_up .* dy_u;
    Vt_up = fulV_up .* dx_v;
    Ut_dw = fulU_dw .* dy_u;
    Vt_dw = fulV_dw .* dx_v;
    Wt_mi = fulW_mi .* dx_u .* dy_v;
    Wt_bo = fulW_bo .* dx_u .* dy_v;
    
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
    
    
    %% SBC Ut
    Ut_up_SBC = zeros(size(Ut_up));
    
    % First
    lat_north_ind_1 = ...
        find(ismember(lat_v, lat_v_SBC_north_repelem(1)));
    lat_south_ind_1 = ...
        find(ismember(lat_v, lat_v_SBC_south_repelem(1)));
    lat_vec_now = lat_north_ind_1:lat_south_ind_1-1;
    lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(1));
    Ut_up_SBC(lat_vec_now, ...
        lon_u_ALLC_repelem_ind) = ...
        Ut_up(lat_vec_now,...
        lon_u_ALLC_repelem_ind);
    % Last
    lat_north_ind_1 = ...
        find(ismember(lat_v, lat_v_SBC_north_repelem(end)));
    lat_south_ind_1 = ...
        find(ismember(lat_v, lat_v_SBC_south_repelem(end)));
    lat_vec_now = lat_north_ind_1:lat_south_ind_1-1;
    lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(end));
    Ut_up_SBC(lat_vec_now, ...
        lon_u_ALLC_repelem_ind) = ...
        Ut_up(lat_vec_now,...
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
        Ut_up_SBC(lat_vec_now, ...
            lon_u_ALLC_repelem_ind) = ...
            Ut_up(lat_vec_now,...
            lon_u_ALLC_repelem_ind);
    end
    
    SBC_Ut_mcps = sum(Ut_up_SBC, 1);
    SBC_Ut = SBC_Ut_mcps .* 10^-6;
    
    
    %% DRC Ut up
    Ut_up_DRC = zeros(size(Ut_up));
    
    % First
    lat_north_ind_1 = ...
        find(ismember(lat_v, lat_v_DRC_north_repelem(1)));
    lat_south_ind_1 = ...
        find(ismember(lat_v, lat_v_DRC_south_repelem(1))) -1;
    lat_vec_now = lat_north_ind_1:lat_south_ind_1;
    lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(1));
    Ut_up_DRC(lat_vec_now, ...
        lon_u_ALLC_repelem_ind) = ...
        Ut_up(lat_vec_now,...
        lon_u_ALLC_repelem_ind);
    % Last
    lat_north_ind_1 = ...
        find(ismember(lat_v, lat_v_DRC_north_repelem(end)));
    lat_south_ind_1 = ...
        find(ismember(lat_v, lat_v_DRC_south_repelem(end))) -1;
    lat_vec_now = lat_north_ind_1:lat_south_ind_1;
    lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(end));
    Ut_up_DRC(lat_vec_now, ...
        lon_u_ALLC_repelem_ind) = ...
        Ut_up(lat_vec_now,...
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
        Ut_up_DRC(lat_vec_now, ...
            lon_u_ALLC_repelem_ind) = ...
            Ut_up(lat_vec_now,...
            lon_u_ALLC_repelem_ind);
    end
    
    DRC_Ut_mcps_up = sum(Ut_up_DRC, 1);
    DRC_Ut_up = DRC_Ut_mcps_up * 10^-6;
    
    
    %% DRC Ut dw
    Ut_dw_DRC = zeros(size(Ut_dw));
    
    % First
    lat_north_ind_1 = ...
        find(ismember(lat_v, lat_v_SBC_north_repelem(1)));
    lat_south_ind_1 = ...
        find(ismember(lat_v, lat_v_DRC_south_repelem(1)));
    lat_vec_now = lat_north_ind_1:lat_south_ind_1-1;
    lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(1));
    Ut_dw_DRC(lat_vec_now, ...
        lon_u_ALLC_repelem_ind) = ...
        Ut_dw(lat_vec_now,...
        lon_u_ALLC_repelem_ind);
    % Last
    lat_north_ind_1 = ...
        find(ismember(lat_v, lat_v_SBC_north_repelem(end)));
    lat_south_ind_1 = ...
        find(ismember(lat_v, lat_v_DRC_south_repelem(end)));
    lat_vec_now = lat_north_ind_1:lat_south_ind_1-1;
    lon_u_ALLC_repelem_ind = find(lon_u == lon_u_ALLC_repelem(end));
    Ut_dw_DRC(lat_vec_now, ...
        lon_u_ALLC_repelem_ind) = ...
        Ut_dw(lat_vec_now,...
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
        Ut_dw_DRC(lat_vec_now, ...
            lon_u_ALLC_repelem_ind) = ...
            Ut_dw(lat_vec_now,...
            lon_u_ALLC_repelem_ind);
    end
    
    DRC_Ut_mcps_dw = sum(Ut_dw_DRC, 1);
    DRC_Ut_dw = DRC_Ut_mcps_dw .* 10^-6;
    
    
    %% DRC Ut
    DRC_Ut = DRC_Ut_up + DRC_Ut_dw;
    
    
    %% SBC Vtn
    Vtn_up_SBC = zeros(1, length(lon_v));
    
    % do the first
    lon_first_idx = 57;
    lat_ind_V_Vtn_up_SBC = ...
        lat_v == lat_v_SBC_north(1);
    V_Vtn_up_SBC_now = ...
        Vt_up(lat_ind_V_Vtn_up_SBC, lon_first_idx);
    V_Vtn_up_SBC = -V_Vtn_up_SBC_now;
    Vtn_up_SBC(lon_first_idx) = V_Vtn_up_SBC;
    
    % and last cases
    lon_last_idx = 312;
    lat_ind_V_Vtn_up_SBC = ...
        lat_v == lat_v_SBC_north(end);
    V_Vtn_up_SBC_now = ...
        Vt_up(lat_ind_V_Vtn_up_SBC, lon_last_idx);
    V_Vtn_up_SBC = -V_Vtn_up_SBC_now;
    Vtn_up_SBC(lon_last_idx) = V_Vtn_up_SBC;
    
    % in between cases
    for jj = 2 : length(lon_u_ALLC)-2
        lon_u_ind = find(lon_u == lon_u_ALLC(jj));
        
        %
        lat_ind_V_Vtn_up_SBC = lat_v == lat_v_SBC_north(jj);
        V_Vtn_up_SBC_now = ...
            Vt_up(lat_ind_V_Vtn_up_SBC, lon_u_ind);
        V_Vtn_up_SBC = -V_Vtn_up_SBC_now;
        
        %
        [U_Vtn_up_SBC_prev, U_Vtn_up_SBC_next] = deal(0);
        
        % if the previous lat is to the south of the current one
        if lat_v_SBC_north(jj-1) < lat_v_SBC_north(jj)
            %
            lat_u_within_ind = find(...
                lat_u < lat_v_SBC_north(jj) & ...
                lat_u > lat_v_SBC_north(jj-1));
            U_Vtn_up_SBC_all_now = ...
                Ut_up(lat_u_within_ind, lon_u_ind);
            U_Vtn_up_SBC_all = U_Vtn_up_SBC_all_now;
            U_Vtn_up_SBC_prev = nansum(U_Vtn_up_SBC_all);
            if isnan(U_Vtn_up_SBC_prev)
                U_Vtn_up_SBC_prev = 0;
            end
        end
        
        % if the next lat is to the south of the current one
        if lat_v_SBC_north(jj+1) < lat_v_SBC_north(jj)
            %
            lat_u_within_ind = find(...
                lat_u < lat_v_SBC_north(jj) & ...
                lat_u > lat_v_SBC_north(jj+1));
            U_Vtn_up_SBC_all_now = ...
                Ut_up(lat_u_within_ind, lon_u_ind+1);
            U_Vtn_up_SBC_all = -U_Vtn_up_SBC_all_now;
            U_Vtn_up_SBC_next = nansum(U_Vtn_up_SBC_all);
        end
        if isnan(U_Vtn_up_SBC_next)
            U_Vtn_up_SBC_next = 0;
        end
        
        %
        Vtn_up_SBC(lon_u_ind) = V_Vtn_up_SBC + ...
            U_Vtn_up_SBC_prev + U_Vtn_up_SBC_next;
    end
    
    
    %% SBC Vts
    Vts_up_SBC = zeros(1, length(lon_v));
    
    % do the first
    lon_first_idx = 57;
    lat_ind_V_Vts_up_SBC = ...
        lat_v == lat_v_SBC_south(1);
    V_Vts_up_SBC_now = ...
        Vt_up(lat_ind_V_Vts_up_SBC, lon_first_idx);
    V_Vts_up_SBC = V_Vts_up_SBC_now;
    Vts_up_SBC(lon_first_idx) = V_Vts_up_SBC;
    
    % and last cases
    lon_last_idx = 312;
    lat_ind_V_Vts_up_SBC = ...
        lat_v == lat_v_SBC_south(end);
    V_Vts_up_SBC_now = ...
        Vt_up(lat_ind_V_Vts_up_SBC, lon_last_idx);
    V_Vts_up_SBC = V_Vts_up_SBC_now;
    Vts_up_SBC(lon_last_idx) = V_Vts_up_SBC;
    
    % in between cases
    for jj = 2 : length(lon_u_ALLC)-2
        lon_u_ind = find(lon_u == lon_u_ALLC(jj));
        
        %
        lat_ind_V_Vts_up_SBC = lat_v == lat_v_SBC_south(jj);
        V_Vts_up_SBC_now = ...
            Vt_up(lat_ind_V_Vts_up_SBC, lon_u_ind);
        V_Vts_up_SBC = V_Vts_up_SBC_now;
        
        %
        [U_Vts_up_SBC_prev, U_Vts_up_SBC_next] = deal(0);
        
        % if the previous lat is to the north of the current one
        if lat_v_SBC_south(jj-1) > lat_v_SBC_south(jj)
            %
            lat_u_within_ind = find(...
                lat_u > lat_v_SBC_south(jj) & ...
                lat_u < lat_v_SBC_south(jj-1));
            U_Vts_up_SBC_all_now = ...
                Ut_up(lat_u_within_ind, lon_u_ind);
            U_Vts_up_SBC_all = U_Vts_up_SBC_all_now;
            U_Vts_up_SBC_prev = nansum(U_Vts_up_SBC_all);
        end
        
        % if the next lat is to the north of the current one
        if lat_v_SBC_south(jj+1) > lat_v_SBC_south(jj)
            %
            lat_u_within_ind = find(...
                lat_u > lat_v_SBC_south(jj) & ...
                lat_u < lat_v_SBC_south(jj+1));
            U_Vts_up_SBC_all_now = ...
                Ut_up(lat_u_within_ind, lon_u_ind+1);
            U_Vts_up_SBC_all = -U_Vts_up_SBC_all_now;
            U_Vts_up_SBC_next = nansum(U_Vts_up_SBC_all);
        end
        
        %
        Vts_up_SBC(lon_u_ind) = ...
            V_Vts_up_SBC + U_Vts_up_SBC_prev + ...
            U_Vts_up_SBC_next;
    end
    
    
    %% SBC Vtc
    SBC_Vtn = Vtn_up_SBC .* 10^-6;
    SBC_Vts = Vts_up_SBC .* 10^-6;
    
    [SBC_Vtnc, SBC_Vtsc] = deal(zeros(1, length(lon_v)));
    SBC_Vtnc(57) = SBC_Vtn(57);
    SBC_Vtsc(57) = SBC_Vts(57);
    
    for jj = 58 : 312
        SBC_Vtnc(jj) = SBC_Vtnc(jj-1) + SBC_Vtn(jj);
        SBC_Vtsc(jj) = SBC_Vtsc(jj-1) + SBC_Vts(jj);
    end
    
    
    %% DRC Vtn up
    Vtn_up_DRC = zeros(1, length(lon_v));
    
    % do the first
    lon_first_idx = 57;
    lat_ind_V_Vtn_up_DRC = ...
        lat_v == lat_v_DRC_north(1);
    V_Vtn_up_DRC_now = ...
        Vt_up(lat_ind_V_Vtn_up_DRC, lon_first_idx);
    V_Vtn_up_DRC = -V_Vtn_up_DRC_now;
    Vtn_up_DRC(lon_first_idx) = V_Vtn_up_DRC;
    
    % and last cases
    lon_last_idx = 312;
    lat_ind_V_Vtn_up_DRC = ...
        lat_v == lat_v_DRC_north(end);
    V_Vtn_up_DRC_now = ...
        Vt_up(lat_ind_V_Vtn_up_DRC, lon_last_idx);
    V_Vtn_up_DRC = -V_Vtn_up_DRC_now;
    Vtn_up_DRC(lon_last_idx) = V_Vtn_up_DRC;
    
    % in between cases
    for jj = length(lon_u_ALLC)-2 : -1 : 2
        lon_u_ind = find(lon_u == lon_u_ALLC(jj));
        
        %
        lat_ind_V_Vtn_up_DRC = lat_v == lat_v_DRC_north(jj);
        V_Vtn_up_DRC_now = ...
            Vt_up(lat_ind_V_Vtn_up_DRC, lon_u_ind);
        V_Vtn_up_DRC = -V_Vtn_up_DRC_now;
        
        %
        [U_Vtn_up_DRC_prev, U_Vtn_up_DRC_next] = deal(0);
        
        % if the previous lat is to the south of the current one
        if lat_v_DRC_north(jj+1) < lat_v_DRC_north(jj)
            %
            lat_u_within_ind = find(...
                lat_u < lat_v_DRC_north(jj) & ...
                lat_u > lat_v_DRC_north(jj+1));
            U_Vtn_up_DRC_all_now = ...
                Ut_up(lat_u_within_ind, lon_u_ind+1);
            U_Vtn_up_DRC_all = -U_Vtn_up_DRC_all_now;
            U_Vtn_up_DRC_prev = nansum(U_Vtn_up_DRC_all);
        end
        
        % if the next lat is to the south of the current one
        if lat_v_DRC_north(jj-1) < lat_v_DRC_north(jj)
            %
            lat_u_within_ind = find(...
                lat_u < lat_v_DRC_north(jj) & ...
                lat_u > lat_v_DRC_north(jj-1));
            U_Vtn_up_DRC_all_now = ...
                Ut_up(lat_u_within_ind, lon_u_ind);
            U_Vtn_up_DRC_all = U_Vtn_up_DRC_all_now;
            U_Vtn_up_DRC_next = nansum(U_Vtn_up_DRC_all);
        end
        
        %
        Vtn_up_DRC(lon_u_ind) = V_Vtn_up_DRC + ...
            U_Vtn_up_DRC_prev + U_Vtn_up_DRC_next;
    end
    
    
    %% DRC Vtn dw
    Vtn_dw_DRC = zeros(1, length(lon_v));
    
    % do the first
    lon_first_idx = 57;
    lat_ind_V_Vtn_dw_DRC = ...
        lat_v == lat_v_SBC_north(1);
    V_Vtn_dw_DRC_now = ...
        Vt_dw(lat_ind_V_Vtn_dw_DRC, lon_first_idx);
    V_Vtn_dw_DRC = -V_Vtn_dw_DRC_now;
    Vtn_dw_DRC(lon_first_idx) = V_Vtn_dw_DRC;
    
    % and last cases
    lon_last_idx = 312;
    lat_ind_V_Vtn_dw_DRC = ...
        lat_v == lat_v_SBC_north(end);
    V_Vtn_dw_DRC_now = ...
        Vt_dw(lat_ind_V_Vtn_dw_DRC, lon_last_idx);
    V_Vtn_dw_DRC = -V_Vtn_dw_DRC_now;
    Vtn_dw_DRC(lon_last_idx) = V_Vtn_dw_DRC;
    
    % in between cases
    for jj = length(lon_u_ALLC)-2 : -1 : 2
        lon_u_ind = find(lon_u == lon_u_ALLC(jj));
        
        %
        lat_ind_V_Vtn_dw_DRC = lat_v == lat_v_SBC_north(jj);
        V_Vtn_dw_DRC_now = ...
            Vt_dw(lat_ind_V_Vtn_dw_DRC, lon_u_ind);
        V_Vtn_dw_DRC = -V_Vtn_dw_DRC_now;
        
        %
        [U_Vtn_dw_DRC_prev, U_Vtn_dw_DRC_next] = deal(0);
        
        % if the previous lat is to the south of the current one
        if lat_v_SBC_north(jj+1) < lat_v_SBC_north(jj)
            %
            lat_u_within_ind = find(...
                lat_u < lat_v_SBC_north(jj) & ...
                lat_u > lat_v_SBC_north(jj+1));
            U_Vtn_dw_DRC_all_now = ...
                Ut_dw(lat_u_within_ind, lon_u_ind+1);
            U_Vtn_dw_DRC_all = -U_Vtn_dw_DRC_all_now;
            U_Vtn_dw_DRC_prev = nansum(U_Vtn_dw_DRC_all);
        end
        
        % if the next lat is to the south of the current one
        if lat_v_SBC_north(jj-1) < lat_v_SBC_north(jj)
            %
            lat_u_within_ind = find(...
                lat_u < lat_v_SBC_north(jj) & ...
                lat_u > lat_v_SBC_north(jj-1));
            U_Vtn_dw_DRC_all_now = ...
                Ut_dw(lat_u_within_ind, lon_u_ind);
            U_Vtn_dw_DRC_all = U_Vtn_dw_DRC_all_now;
            U_Vtn_dw_DRC_next = nansum(U_Vtn_dw_DRC_all);
        end
        
        %
        Vtn_dw_DRC(lon_u_ind) = V_Vtn_dw_DRC + ...
            U_Vtn_dw_DRC_prev + U_Vtn_dw_DRC_next;
    end
    
    for n = 1 : length(lon_u)-1
        if Vtn_dw_DRC(n) ~= 0
            error('Vtn_dw_DRC ~= 0 somewhere :( !')
        end
    end
    disp('Vtn_dw_DRC = 0 everywhere :)')
    
    
    %% DRC Vts up
    Vts_up_DRC = zeros(1, length(lon_v));
    
    % do the first
    lon_first_idx = 57;
    lat_ind_V_Vts_up_DRC = ...
        lat_v == lat_v_DRC_south(1);
    V_Vts_up_DRC_now = ...
        Vt_up(lat_ind_V_Vts_up_DRC, lon_first_idx);
    V_Vts_up_DRC = V_Vts_up_DRC_now;
    Vts_up_DRC(lon_first_idx) = V_Vts_up_DRC;
    
    % and last cases
    lon_last_idx = 312;
    lat_ind_V_Vts_up_DRC = ...
        lat_v == lat_v_DRC_south(end);
    V_Vts_up_DRC_now = ...
        Vt_up(lat_ind_V_Vts_up_DRC, lon_last_idx);
    V_Vts_up_DRC = V_Vts_up_DRC_now;
    Vts_up_DRC(lon_last_idx) = V_Vts_up_DRC;
    
    % in between cases
    for jj = length(lon_u_ALLC)-2 : -1 : 2
        lon_u_ind = find(lon_u == lon_u_ALLC(jj));
        
        %
        lat_ind_V_Vts_up_DRC = lat_v == lat_v_DRC_south(jj);
        V_Vts_up_DRC_now = ...
            Vt_up(lat_ind_V_Vts_up_DRC, lon_u_ind);
        V_Vts_up_DRC = V_Vts_up_DRC_now;
        
        %
        [U_Vts_up_DRC_prev, U_Vts_up_DRC_next] = deal(0);
        
        % if the previous lat is to the north of the current one
        if lat_v_DRC_south(jj+1) > lat_v_DRC_south(jj)
            %
            lat_u_within_ind = find(...
                lat_u > lat_v_DRC_south(jj) & ...
                lat_u < lat_v_DRC_south(jj+1));
            U_Vts_up_DRC_all_now = ...
                Ut_up(lat_u_within_ind, lon_u_ind+1);
            U_Vts_up_DRC_all = -U_Vts_up_DRC_all_now;
            U_Vts_up_DRC_prev = nansum(U_Vts_up_DRC_all);
        end
        
        % if the next lat is to the north of the current one
        if lat_v_DRC_south(jj-1) > lat_v_DRC_south(jj)
            %
            lat_u_within_ind = find(...
                lat_u > lat_v_DRC_south(jj) & ...
                lat_u < lat_v_DRC_south(jj-1));
            U_Vts_up_DRC_all_now = ...
                Ut_up(lat_u_within_ind, lon_u_ind);
            U_Vts_up_DRC_all = U_Vts_up_DRC_all_now;
            U_Vts_up_DRC_next = nansum(U_Vts_up_DRC_all);
        end
        
        %
        Vts_up_DRC(lon_u_ind) = V_Vts_up_DRC + ...
            U_Vts_up_DRC_prev + U_Vts_up_DRC_next;
    end
    
    
    %% DRC Vts dw
    Vts_dw_DRC = zeros(1, length(lon_v));
    
    % do the first
    lon_first_idx = 57;
    lat_ind_V_Vts_dw_DRC = ...
        lat_v == lat_v_DRC_south(1);
    V_Vts_dw_DRC_now = ...
        Vt_dw(lat_ind_V_Vts_dw_DRC, lon_first_idx);
    V_Vts_dw_DRC = V_Vts_dw_DRC_now;
    Vts_dw_DRC(lon_first_idx) = V_Vts_dw_DRC;
    
    % and last cases
    lon_last_idx = 312;
    lat_ind_V_Vts_dw_DRC = ...
        lat_v == lat_v_DRC_south(end);
    V_Vts_dw_DRC_now = ...
        Vt_dw(lat_ind_V_Vts_dw_DRC, lon_last_idx);
    V_Vts_dw_DRC = V_Vts_dw_DRC_now;
    Vts_dw_DRC(lon_last_idx) = V_Vts_dw_DRC;
    
    % in between cases
    for jj = length(lon_u_ALLC)-2 : -1 : 2
        lon_u_ind = find(lon_u == lon_u_ALLC(jj));
        
        %
        lat_ind_V_Vts_dw_DRC = lat_v == lat_v_DRC_south(jj);
        V_Vts_dw_DRC_now = ...
            Vt_dw(lat_ind_V_Vts_dw_DRC, lon_u_ind);
        V_Vts_dw_DRC = V_Vts_dw_DRC_now;
        
        %
        [U_Vts_dw_DRC_prev, U_Vts_dw_DRC_next] = deal(0);
        
        % if the previous lat is to the north of the current one
        if lat_v_DRC_south(jj+1) > lat_v_DRC_south(jj)
            %
            lat_u_within_ind = find(...
                lat_u > lat_v_DRC_south(jj) & ...
                lat_u < lat_v_DRC_south(jj+1));
            U_Vts_dw_DRC_all_now = ...
                Ut_dw(lat_u_within_ind, lon_u_ind+1);
            U_Vts_dw_DRC_all = -U_Vts_dw_DRC_all_now;
            U_Vts_dw_DRC_prev = nansum(U_Vts_dw_DRC_all);
        end
        
        % if the next lat is to the north of the current one
        if lat_v_DRC_south(jj-1) > lat_v_DRC_south(jj)
            %
            lat_u_within_ind = find(...
                lat_u > lat_v_DRC_south(jj) & ...
                lat_u < lat_v_DRC_south(jj-1));
            U_Vts_dw_DRC_all_now = ...
                Ut_dw(lat_u_within_ind, lon_u_ind);
            U_Vts_dw_DRC_all = U_Vts_dw_DRC_all_now;
            U_Vts_dw_DRC_next = nansum(U_Vts_dw_DRC_all);
        end
        
        %
        Vts_dw_DRC(lon_u_ind) = ...
            V_Vts_dw_DRC + ...
            U_Vts_dw_DRC_prev + U_Vts_dw_DRC_next;
    end
    
    
    %% DRC Vtc
    DRC_Vtn = Vtn_up_DRC .* 10^-6;
    DRC_Vts_up = Vts_up_DRC .* 10^-6;
    DRC_Vts_dw = Vts_dw_DRC .* 10^-6;
    
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
    % Get Ut_up_SBC Vt_up_SBC for W
    Wt_mi_SBC_for_W = NaN(size(fulV_up));
    
    nc = 0;
    for n = 57 : 312
        nc = nc + 1;
        lat_v_SBC_north_ind = find(lat_v==lat_v_SBC_north(nc));
        lat_v_SBC_south_ind = find(lat_v==lat_v_SBC_south(nc));
        lat_v_SBC_ind = lat_v_SBC_north_ind : lat_v_SBC_south_ind;
        
        Wt_mi_SBC_for_W(lat_v_SBC_ind,n) = Wt_mi(lat_v_SBC_ind,n);
    end
    
    SBC_Wt_mcps_real = nansum(Wt_mi_SBC_for_W,1);
    SBC_Wt = SBC_Wt_mcps_real * 10^-6;
    
    
    %% DRC Wt
    DRC_Wt = -SBC_Wt;
    
    
    %% BOT Wt
    % Get Ut_up_SBC Vt_up_SBC for W
    Wt_bo_BOT_for_W = NaN(size(fulV_up));
    
    nc = 0;
    for n = 57 : 312
        nc = nc + 1;
        lat_v_BOT_north_ind = find(lat_v==lat_v_SBC_north(nc));
        lat_v_BOT_south_ind = find(lat_v==lat_v_DRC_south(nc));
        lat_v_BOT_ind = lat_v_BOT_north_ind : lat_v_BOT_south_ind;
        
        Wt_bo_BOT_for_W(lat_v_BOT_ind,n) = Wt_bo(lat_v_BOT_ind,n);
    end
    
    BOT_Wt_mcps_real = nansum(Wt_bo_BOT_for_W,1);
    BOT_Wt = BOT_Wt_mcps_real * 10^-6;
    
    
    %% Wtc
    [SBC_Wtc, DRC_Wtc, BOT_Wtc] = deal(zeros(1, length(lon_v)));
    
    SBC_Wtc(57) = SBC_Wt(57);
    for jj = 58 : 312
        SBC_Wtc(jj) = SBC_Wtc(jj-1) + SBC_Wt(jj);
    end
    
    DRC_Wtc(312) = DRC_Wt(312);
    BOT_Wtc(312) = BOT_Wt(312);
    for jj = 311 : -1 : 57
        DRC_Wtc(jj) = DRC_Wtc(jj+1) + DRC_Wt(jj);
        BOT_Wtc(jj) = BOT_Wtc(jj+1) + BOT_Wt(jj);
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
            -(DRC_Vtnc(jj) + DRC_Vtsc(jj) + DRC_Wtc(jj) + BOT_Wtc(jj));
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
    plot(lon_u, SBC_Ut, 'k-', 'linewidth', 1)
    hold all
    plot(lon_u, SBC_Ut_star, 'k--', 'linewidth', 1)
    plot(lon_v, SBC_Vtnc, 'g-', 'linewidth', 1)
    plot(lon_v, SBC_Vtsc, 'r-', 'linewidth', 1)
    plot(lon_v, SBC_Wtc, 'b-', 'linewidth', 1)
    legend(...
        'SBC $U_{t}$', ...
        'SBC $U_{t}^{*}$', ...
        'SBC $V_{t}nc$', ...
        'SBC $V_{t}sc$', ...
        'SBC $W_{t}c$')
    axis([115 147 -3 3])
    title([MTH{t} ... 
        ' KDS75 SBC transports ($Sv$) VS longitude'])
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
    plot(lon_u, DRC_Ut, 'k-', 'linewidth', 1)
    hold all
    plot(lon_u, DRC_Ut_star, 'k--', 'linewidth', 1)
    plot(lon_v, DRC_Vtnc, 'r-', 'linewidth', 1)
    plot(lon_v, DRC_Vtsc, 'm-', 'linewidth', 1)
    plot(lon_v, DRC_Vtsc_up, 'm:', 'linewidth', 1)
    plot(lon_v, DRC_Vtsc_dw, 'm--', 'linewidth', 1)
    plot(lon_v, DRC_Wtc, 'b-', 'linewidth', 1)
    plot(lon_v, BOT_Wtc, 'y-', 'linewidth', 1)
    legend(...
        'DRC $U_{t}$', ...
        'DRC $U_{t}^{*}$', ...
        'DRC $V_{t}nc$', ...
        'DRC $V_{t}sc$', ...
        'DRC $V_{t}sc$ up', ...
        'DRC $V_{t}sc$ dw', ...
        'DRC $W_{t}c$', ...
        'BOT $W_{t}c$')
    axis([115 147 -17 10])
    title([MTH{t} ... 
        ' KDS75 DRC transports ($Sv$) VS longitude'])
    grid
    set(gca,'layer','top','color',fig_color,...
        'fontsize',font_size,'tickdir','out')
    
    outputls = ls(figures_path);
    scriptname = mfilename;
    if ~contains(outputls, scriptname)
        mkdir(figures_path, scriptname)
    end
    export_fig(fig, ...
        [figures_path mfilename '/' scriptname(1:3) ...
        '_fig' num2str(fig_n) '_' num2str(t)], ...
        '-m1')
    
    
    %% save all that stuff !
    KDau_fcrt.SBC_Ut.(MTH{t}) = SBC_Ut;
    KDau_fcrt.SBC_Ut_star.(MTH{t}) = SBC_Ut_star;
    KDau_fcrt.SBC_Vtnc.(MTH{t}) = SBC_Vtnc;
    KDau_fcrt.SBC_Vtsc.(MTH{t}) = SBC_Vtsc;
    KDau_fcrt.SBC_Wtc.(MTH{t}) = SBC_Wtc;
    KDau_fcrt.DRC_Ut.(MTH{t}) = DRC_Ut;
    KDau_fcrt.DRC_Ut_star.(MTH{t}) = DRC_Ut_star;
    KDau_fcrt.DRC_Vtnc.(MTH{t}) = DRC_Vtnc;
    KDau_fcrt.DRC_Vtsc.(MTH{t}) = DRC_Vtsc;
    KDau_fcrt.DRC_Vtsc_up.(MTH{t}) = DRC_Vtsc_up;
    KDau_fcrt.DRC_Vtsc_dw.(MTH{t}) = DRC_Vtsc_dw;
    KDau_fcrt.DRC_Wtc.(MTH{t}) = DRC_Wtc;
    KDau_fcrt.BOT_Wtc.(MTH{t}) = BOT_Wtc;
    disp([MTH{t} ' OK!'])
end


%%
save([data_path 'SACS_data/KDau_fcrt'], 'KDau_fcrt')
disp(['KDau_fcrt DONE'])

