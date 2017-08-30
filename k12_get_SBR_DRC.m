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
lat_v_SBC_north = aus8_currents.lat_v_SBC_north;
lat_v_SBC_south = aus8_currents.lat_v_SBC_south;
lat_v_DRC_north = aus8_currents.lat_v_DRC_north;
lat_v_DRC_south = aus8_currents.lat_v_DRC_south;
lon_u_ALLC = aus8_currents.lon_u_ALLC;
Months = aus8_coor.Months;
Months{13} = 'mean';


%% Prep
for t = 1 : 13
    %
    U_prime_up = KDau_currents.ztop_to_zmid.U_prime.(Months{t});
    V_prime_up = KDau_currents.ztop_to_zmid.V_prime.(Months{t});
    U_g_prime_dw = KDau_currents.zmid_to_zbot.U_g_prime.(Months{t});
    V_g_prime_dw = KDau_currents.zmid_to_zbot.V_g_prime.(Months{t});
    
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
    KDau_U_prime.(Months{t})(isnan(KDau_U_prime.(Months{t}))) = 0;
    KDau_V_prime.(Months{t})(isnan(KDau_V_prime.(Months{t}))) = 0;
    du = KDau_U_prime.(Months{t})(:,2:end) - ...
        KDau_U_prime.(Months{t})(:,1:end-1);
    dudx = du ./ dx_u;
    lat_v_repmat = repmat(lat_v,1,length(lon_v));
    dv = KDau_V_prime.(Months{t})(1:end-1,:).* ...
        cos(lat_v_repmat(1:end-1,:) * pi180) - ...
        KDau_V_prime.(Months{t})(2:end,:).* ...
        cos(lat_v_repmat(2:end,:) * pi180);
    lat_u_repmat = repmat(lat_u,1,length(lon_v));
    dvdy = 1./cos(lat_u_repmat * pi180).*dv./dy_v;
    div_UV_prime.(Months{t}) = dudx + dvdy;
    
    
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
    DRC_Ut_up = DRC_Ut_mcps_up * 10^-6;
    
    
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
                Ut_prime_up(lat_u_within_ind, lon_u_ind+1);
            U_Vtn_prime_up_DRC_all = -U_Vtn_prime_up_DRC_all_now;
            U_Vtn_prime_up_DRC_prev = nansum(U_Vtn_prime_up_DRC_all);
        end
        
        % if the next lat is to the south of the current one
        if lat_v_DRC_north(jj-1) < lat_v_DRC_north(jj)
            %
            lat_u_within_ind = find(...
                lat_u < lat_v_DRC_north(jj) & ...
                lat_u > lat_v_DRC_north(jj-1));
            U_Vtn_prime_up_DRC_all_now = ...
                Ut_prime_up(lat_u_within_ind, lon_u_ind);
            U_Vtn_prime_up_DRC_all = U_Vtn_prime_up_DRC_all_now;
            U_Vtn_prime_up_DRC_next = nansum(U_Vtn_prime_up_DRC_all);
        end
        
        %
        Vtn_prime_up_DRC(lon_u_ind) = V_Vtn_prime_up_DRC + ...
            U_Vtn_prime_up_DRC_prev + U_Vtn_prime_up_DRC_next;
    end
    
    
    %% DRC Vtn dw
    Vtn_g_prime_dw_DRC = zeros(1, length(lon_v));
    
    % do the first
    lon_first_idx = 57;
    lat_ind_V_Vtn_g_prime_dw_DRC = ...
        lat_v == lat_v_SBC_north(1);
    V_Vtn_g_prime_dw_DRC_now = ...
        Vt_g_prime_dw(lat_ind_V_Vtn_g_prime_dw_DRC, lon_first_idx);
    V_Vtn_g_prime_dw_DRC = -V_Vtn_g_prime_dw_DRC_now;
    Vtn_g_prime_dw_DRC(lon_first_idx) = V_Vtn_g_prime_dw_DRC;
    
    % and last cases
    lon_last_idx = 312;
    lat_ind_V_Vtn_g_prime_dw_DRC = ...
        lat_v == lat_v_SBC_north(end);
    V_Vtn_g_prime_dw_DRC_now = ...
        Vt_g_prime_dw(lat_ind_V_Vtn_g_prime_dw_DRC, lon_last_idx);
    V_Vtn_g_prime_dw_DRC = -V_Vtn_g_prime_dw_DRC_now;
    Vtn_g_prime_dw_DRC(lon_last_idx) = V_Vtn_g_prime_dw_DRC;
    
    % in between cases
    for jj = length(lon_u_ALLC)-2 : -1 : 2
        lon_u_ind = find(lon_u == lon_u_ALLC(jj));
        
        %
        lat_ind_V_Vtn_g_prime_dw_DRC = lat_v == lat_v_SBC_north(jj);
        V_Vtn_g_prime_dw_DRC_now = ...
            Vt_g_prime_dw(lat_ind_V_Vtn_g_prime_dw_DRC, lon_u_ind);
        V_Vtn_g_prime_dw_DRC = -V_Vtn_g_prime_dw_DRC_now;
        
        %
        [U_Vtn_g_prime_dw_DRC_prev, U_Vtn_g_prime_dw_DRC_next] = deal(0);
        
        % if the previous lat is to the south of the current one
        if lat_v_SBC_north(jj+1) < lat_v_SBC_north(jj)
            %
            lat_u_within_ind = find(...
                lat_u < lat_v_SBC_north(jj) & ...
                lat_u > lat_v_SBC_north(jj+1));
            U_Vtn_g_prime_dw_DRC_all_now = ...
                Ut_g_prime_dw(lat_u_within_ind, lon_u_ind+1);
            U_Vtn_g_prime_dw_DRC_all = -U_Vtn_g_prime_dw_DRC_all_now;
            U_Vtn_g_prime_dw_DRC_prev = nansum(U_Vtn_g_prime_dw_DRC_all);
        end
        
        % if the next lat is to the south of the current one
        if lat_v_SBC_north(jj-1) < lat_v_SBC_north(jj)
            %
            lat_u_within_ind = find(...
                lat_u < lat_v_SBC_north(jj) & ...
                lat_u > lat_v_SBC_north(jj-1));
            U_Vtn_g_prime_dw_DRC_all_now = ...
                Ut_g_prime_dw(lat_u_within_ind, lon_u_ind);
            U_Vtn_g_prime_dw_DRC_all = U_Vtn_g_prime_dw_DRC_all_now;
            U_Vtn_g_prime_dw_DRC_next = nansum(U_Vtn_g_prime_dw_DRC_all);
        end
        
        %
        Vtn_g_prime_dw_DRC(lon_u_ind) = V_Vtn_g_prime_dw_DRC + ...
            U_Vtn_g_prime_dw_DRC_prev + U_Vtn_g_prime_dw_DRC_next;
    end
    
    for n = 1 : length(lon_u)-1
        if Vtn_g_prime_dw_DRC(n) ~= 0
            error('Vtn_g_prime_dw_DRC = 0 somewhere :( !')
        end
    end
    disp('Vtn_g_prime_dw_DRC = 0 everywhere :)')
    
    
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
        Vts_prime_up_DRC(lon_u_ind) = V_Vts_prime_up_DRC + ...
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
    w_prime_up_SBC_real = div_UV_prime_up_SBC;
    w_prime_up_SBC_leak = div_UV_prime_up_SBC;
    
    for ii = 1 : length(lat_u)
        for jj = 1 : length(lon_v)
            if isfinite(div_UV_prime_up_SBC(ii,jj))
                if div_UV_prime_up_SBC(ii,jj) == ...
                        div_UV_prime.(Months{t})(ii,jj)
                    w_prime_up_SBC_real(ii,jj) = 0;
                    
                else
                    w_prime_up_SBC_leak(ii,jj) = 0;
                end
            end
        end
    end
    
    W_prime_up_mcps_SBC_real = w_prime_up_SBC_real .* dy_v .* dx_u;
    SBC_Wt_mcps_real = nansum(W_prime_up_mcps_SBC_real,1);
    SBC_Wt_real = SBC_Wt_mcps_real * 10^-6;
    
    W_prime_up_mcps_SBC_leak = w_prime_up_SBC_leak .* dy_v .* dx_u;
    SBC_Wt_mcps_leak = nansum(W_prime_up_mcps_SBC_leak,1);
    SBC_Wt_leak = SBC_Wt_mcps_leak * 10^-6;
    
    SBC_Wt = SBC_Wt_real + SBC_Wt_leak;
    
    
    %% DRC Wt
    % Get Ut_prime_up_SBC Vt_prime_up_SBC for W
    V_g_prime_dw_DRC_for_W = NaN(size(V_prime_up));
    U_g_prime_dw_DRC_for_W = NaN(size(U_prime_up));
    
    nc = 257;
    for n = 313 : -1 : 58
        nc = nc - 1;
        lat_v_DRC_north_ind = find(lat_v==lat_v_SBC_north(nc));
        lat_v_DRC_south_ind = find(lat_v==lat_v_DRC_south(nc));
        lat_v_DRC_ind = lat_v_DRC_north_ind : lat_v_DRC_south_ind;
        
        lat_u_DRC_north_ind = find(lat_u<lat_v_SBC_north(nc), 1, 'first');
        lat_u_DRC_south_ind = find(lat_u>lat_v_DRC_south(nc), 1, 'last');
        lat_u_DRC_ind = lat_u_DRC_north_ind : lat_u_DRC_south_ind;
        
        V_g_prime_dw_DRC_for_W(lat_v_DRC_ind,n-1) = ...
            V_g_prime_dw(lat_v_DRC_ind,n-1);
        U_g_prime_dw_DRC_for_W(lat_u_DRC_ind,n-1:n) = ...
            U_g_prime_dw(lat_u_DRC_ind,n-1:n);
    end
    
    du = U_g_prime_dw_DRC_for_W(:,2:end) - ...
        U_g_prime_dw_DRC_for_W(:,1:end-1);
    dudx = du ./ dx_u;
    dv = V_g_prime_dw_DRC_for_W(1:end-1,:).* ...
        cos(lat_v_repmat(1:end-1,:) * pi180) - ...
        V_g_prime_dw_DRC_for_W(2:end,:).* ...
        cos(lat_v_repmat(2:end,:) * pi180);
    dvdy = 1./cos(lat_u_repmat * pi180).*dv./dy_v;
    div_UV_g_prime_dw_DRC = dudx + dvdy;
    
    w_g_prime_dw_DRC_real = -w_prime_up_SBC_real;
    
    %%% calculate UV prime in the whole upper layer
    du = U_prime_up(:,2:end) - ...
        U_prime_up(:,1:end-1);
    dudx = du ./ dx_u;
    dv = V_prime_up(1:end-1,:).* ...
        cos(lat_v_repmat(1:end-1,:) * pi180) - ...
        V_prime_up(2:end,:).* ...
        cos(lat_v_repmat(2:end,:) * pi180);
    dvdy = 1./cos(lat_u_repmat * pi180).*dv./dy_v;
    div_UV_prime_up_all = dudx + dvdy;
    %%%
    
    div_UV_for_this = div_UV_g_prime_dw_DRC + div_UV_prime_up_all;
    div_UV_for_this(w_g_prime_dw_DRC_real==0) = NaN;
    w_prime_dw_DRC_leak = div_UV_for_this;
    
    W_g_prime_dw_mcps_DRC_real = w_g_prime_dw_DRC_real .* dy_v .* dx_u;
    DRC_Wt_mcps_real = nansum(W_g_prime_dw_mcps_DRC_real,1);
    DRC_Wt_real = DRC_Wt_mcps_real * 10^-6;
    
    W_g_prime_dw_mcps_DRC_leak = w_prime_dw_DRC_leak .* dy_v .* dx_u;
    DRC_Wt_mcps_leak = nansum(W_g_prime_dw_mcps_DRC_leak,1);
    DRC_Wt_leak = DRC_Wt_mcps_leak * 10^-6;
    
    DRC_Wt = DRC_Wt_real + DRC_Wt_leak;
    
    
    %% Wtc
    [SBC_Wtc_real, SBC_Wtc_leak, DRC_Wtc_real, DRC_Wtc_leak] = ...
        deal(zeros(1, length(lon_v)));
    
    SBC_Wtc_real(57) = SBC_Wt_real(57);
    SBC_Wtc_leak(57) = SBC_Wt_leak(57);
    for jj = 58 : 312
        SBC_Wtc_real(jj) = SBC_Wtc_real(jj-1) + SBC_Wt_real(jj);
        SBC_Wtc_leak(jj) = SBC_Wtc_leak(jj-1) + SBC_Wt_leak(jj);
    end
    SBC_Wtc = SBC_Wtc_real + SBC_Wtc_leak;
    
    DRC_Wtc_real(312) = DRC_Wt_real(312);
    DRC_Wtc_leak(312) = DRC_Wt_leak(312);
    for jj = 311 : -1 : 57
        DRC_Wtc_real(jj) = DRC_Wtc_real(jj+1) + DRC_Wt_real(jj);
        DRC_Wtc_leak(jj) = DRC_Wtc_leak(jj+1) + DRC_Wt_leak(jj);
    end
    DRC_Wtc = DRC_Wtc_real + DRC_Wtc_leak;
    
    
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
    plot(lon_u, SBC_Ut, 'k-', 'linewidth', 1)
    hold all
    plot(lon_u, SBC_Ut_star, 'k--', 'linewidth', 1)
    plot(lon_v, SBC_Vtnc, 'g-', 'linewidth', 1)
    plot(lon_v, SBC_Vtsc, 'r-', 'linewidth', 1)
    plot(lon_v, SBC_Wtc_real, 'b-', 'linewidth', 1)
    plot(lon_v, SBC_Wtc_leak, 'b:', 'linewidth', 1)
    legend(...
        'SBC $U_{t}$', ...
        'SBC $U_{t}^{*}$', ...
        'SBC $V_{t}nc$', ...
        'SBC $V_{t}sc$', ...
        'SBC $W_{t}c$ real', ...
        'SBC $W_{t}c$ leak')
    axis([115 147 -3 3])
    title([Months{t} ... 
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
    plot(lon_v, DRC_Wtc_real, 'b-', 'linewidth', 1)
    plot(lon_v, DRC_Wtc_leak, 'b:', 'linewidth', 1)
    legend(...
        'DRC $U_{t}$', ...
        'DRC $U_{t}^{*}$', ...
        'DRC $V_{t}nc$', ...
        'DRC $V_{t}sc$', ...
        'DRC $V_{t}sc$ up', ...
        'DRC $V_{t}sc$ dw', ...
        'DRC $W_{t}c$ real', ...
        'DRC $W_{t}c$ leak')
    axis([115 147 -17 17])
    title([Months{t} ... 
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
    
    
    %%
    close all
    fig_n = 2;
    rowcols = [2 2];
    rowcols_size = [16 6]; % cm
    margs = [1.5 1 1 1]; % cm
    gaps = [1.5 2]; % cm
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
    plot(lon_u, SBC_Ut-SBC_Ut_star, 'k-', 'linewidth', 1)
    title([Months{t} ...
        ' KDS75 SBC $U_{t}-U_{t}^{*}$ ($Sv$) leak included VS longitude'])
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
    plot(lon_u, DRC_Ut-DRC_Ut_star, 'k-', 'linewidth', 1)
    title([Months{t} ...
        ' KDS75 DRC $U_{t}-U_{t}^{*}$ ($Sv$) leak included VS longitude'])
    grid
    set(gca,'layer','top','color',fig_color,...
        'fontsize',font_size,'tickdir','out')
    
    sp = 3;
    axes('Units','centimeters', ...
        'Position',[...
        (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)), ...
        (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
        x_sp, ...
        y_sp])
    plot(lon_v, SBC_Wtc_leak, 'b:', 'linewidth', 1)
    title([Months{t} ...
        ' KDS75 SBC $W_{t}c$ leak ($Sv$) VS longitude'])
    grid
    set(gca,'layer','top','color',fig_color,...
        'fontsize',font_size,'tickdir','out')
    
    sp = 4;
    axes('Units','centimeters', ...
        'Position',[...
        (marg_l+x_sp*(cm(sp)-1)+gap_w*(cm(sp)-1)), ...
        (marg_b+y_sp*(rm(sp)-1)+gap_h*(rm(sp)-1)), ...
        x_sp, ...
        y_sp])
    plot(lon_v, DRC_Wtc_leak, 'b:', 'linewidth', 1)
    title('KDS75 DRC $W_{t}c$ leak ($Sv$) VS longitude')
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
    close
    
    
    %% save all that stuff !
    KDau_currents.SBC_Ut.(Months{t}) = SBC_Ut;
    KDau_currents.SBC_Ut_star.(Months{t}) = SBC_Ut_star;
    KDau_currents.SBC_Vtnc.(Months{t}) = SBC_Vtnc;
    KDau_currents.SBC_Vtsc.(Months{t}) = SBC_Vtsc;
    KDau_currents.SBC_Wtc_real.(Months{t}) = SBC_Wtc_real;
    KDau_currents.SBC_Wtc_leak.(Months{t}) = SBC_Wtc_leak;
    KDau_currents.DRC_Ut.(Months{t}) = DRC_Ut;
    KDau_currents.DRC_Ut_star.(Months{t}) = DRC_Ut_star;
    KDau_currents.DRC_Vtnc.(Months{t}) = DRC_Vtnc;
    KDau_currents.DRC_Vtsc.(Months{t}) = DRC_Vtsc;
    KDau_currents.DRC_Vtsc_up.(Months{t}) = DRC_Vtsc_up;
    KDau_currents.DRC_Vtsc_dw.(Months{t}) = DRC_Vtsc_dw;
    KDau_currents.DRC_Wtc_real.(Months{t}) = DRC_Wtc_real;
    KDau_currents.DRC_Wtc_leak.(Months{t}) = DRC_Wtc_leak;
    disp([Months{t} ' OK!'])
end


%%
save([data_path 'SACS_data/KDau_currents'], 'KDau_currents')
disp(['KDau_currents DONE'])

