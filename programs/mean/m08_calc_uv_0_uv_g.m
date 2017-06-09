%% calculate geostrophic velocity assuming f is constant
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_dynh_0'])

lat = aus8_coor.lat;
lon = aus8_coor.lon;
depth = aus8_coor.depth;
Months = aus8_coor.Months;

a = 6371000; % Earth's radius in meters
pi180 = pi/180;
lat_v = lat;
lat_u = (lat(1:end-1) + lat(2:end))/2;
lon_u = lon;
lon_v = (lon(1:end-1) + lon(2:end))/2;

aus8_coor.a = a;
aus8_coor.pi180 = pi180;
aus8_coor.lat_v = lat_v;
aus8_coor.lat_u = lat_u;
aus8_coor.lon_u = lon_u;
aus8_coor.lon_v = lon_v;
save([data_path 'SACS_data/aus8_coor'], 'aus8_coor')


%% u_0
aus8_u_0.mean = NaN(length(lat_u),length(lon_u),length(depth));
for t = 1 : 12
    aus8_u_0.(Months{t}) = ...
        NaN(length(lat_u),length(lon_u),length(depth));
end
f_u = -gsw_f(lat_u);
f_u_repmat = repmat(f_u, 1, length(depth));
for n = 1 : length(lon_u)
    % dynamic height at that longitude
    dynh_now = squeeze(aus8_dynh_0.mean(:,n,:));
    
    % dynamic height difference
    dynh_diff_now = dynh_now(1:end-1,:) - dynh_now(2:end,:);
    
    % latitude distance
    dy_now = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
    dy_now_repmat = repmat(dy_now,1,length(depth));
    
    % u derivation
    u_now = dynh_diff_now ./ (dy_now_repmat .* f_u_repmat);
    
    % save the current u and turn into cm/s
    aus8_u_0.mean(:,n,:) = u_now;
    
    % monthly
    for t = 1 : 12
        dynh_now = squeeze(aus8_dynh_0.(Months{t})(:,n,:));
        dynh_diff_now = dynh_now(1:end-1,:) - dynh_now(2:end,:);
        u_now = dynh_diff_now ./ (dy_now_repmat .* f_u_repmat);
        aus8_u_0.(Months{t})(:,n,:) = u_now;
    end
    fprintf('lon_u = %7.3f \n', lon_u(n))
end


% v_0
aus8_v_0.mean = NaN(length(lat_v),length(lon_v),length(depth));
Months = {...
    'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
for t = 1 : 12
    aus8_v_0.(Months{t}) = ...
        NaN(length(lat_v),length(lon_v),length(depth));
end
for m = 1 : length(lat_v)
    f_v = gsw_f(lat_v(m));
    f_v_repmat = repmat(f_v, length(lon_v), length(depth));
    
    % dynamic height at that latitude
    dynh_now = squeeze(aus8_dynh_0.mean(m,:,:));
    
    % dynamic height difference
    dynh_diff_now = ...
        dynh_now(2:end,:) - dynh_now(1:end-1,:);
    
    % longitude distance
    dx_now = a * cos(lat_v(m) * pi180) * ...
        (lon_u(2:end) - lon_u(1:end-1)) * pi180;
    dx_now_repmat = repmat(dx_now',1,length(depth));
    
    % v derivation
    v_now = dynh_diff_now ./ (dx_now_repmat .* f_v_repmat);
    
    % save the current v and convert into cm/s
    aus8_v_0.mean(m,:,:) = v_now;
    
    % monthly
    for t = 1 : 12
        dynh_now = squeeze(aus8_dynh_0.(Months{t})(m,:,:));
        dynh_diff_now = ...
            dynh_now(2:end,:) - dynh_now(1:end-1,:);
        v_now = dynh_diff_now ./ (dx_now_repmat .* f_v_repmat);
        aus8_v_0.(Months{t})(m,:,:) = v_now;
    end
    fprintf('lat_v = %7.3f \n', lat(m))
end


%% Adjust the geostrophic velocities
aus8_u_0_botm.mean = NaN(size(length(lat_u),length(lon_u)));
aus8_v_0_botm.mean = NaN(size(length(lat_v),length(lon_v)));
for t = 1 : 12
    aus8_u_0_botm.(Months{t}) = NaN(length(lat_u),length(lon_u));
    aus8_v_0_botm.(Months{t}) = NaN(length(lat_v),length(lon_v));
end

for m = 1 : length(lat_u)
    for n = 1 : length(lon_u)
        last_finite = find(isfinite(aus8_u_0.mean(m,n,:)),1,'last');
        if ~isempty(last_finite)
            aus8_u_0_botm.mean(m,n) = aus8_u_0.mean(m,n,last_finite);
            for t = 1 : 12
                aus8_u_0_botm.(Months{t})(m,n) = ...
                    aus8_u_0.(Months{t})(m,n,last_finite);
            end
        end
    end
end
aus8_u_g_dum.mean = aus8_u_0.mean - ...
    repmat(aus8_u_0_botm.mean, [1 1 length(depth)]);
for t = 1 : 12
    aus8_u_g_dum.(Months{t}) = aus8_u_0.(Months{t}) - ...
        repmat(aus8_u_0_botm.(Months{t}), [1 1 length(depth)]);
end

for m = 1 : length(lat_v)
    for n = 1 : length(lon_v)
        last_finite = find(isfinite(aus8_v_0.mean(m,n,:)),1,'last');
        if ~isempty(last_finite)
            aus8_v_0_botm.mean(m,n) = aus8_v_0.mean(m,n,last_finite);
            for t = 1 : 12
                aus8_v_0_botm.(Months{t})(m,n) = ...
                    aus8_v_0.(Months{t})(m,n,last_finite);
            end
        end
    end
end
aus8_v_g_dum.mean = aus8_v_0.mean - ...
    repmat(aus8_v_0_botm.mean, [1 1 length(depth)]);
for t = 1 : 12
    aus8_v_g_dum.(Months{t}) = aus8_v_0.(Months{t}) - ...
        repmat(aus8_v_0_botm.(Months{t}), [1 1 length(depth)]);
end


%% get the uv_g velocities at mid depth point
depth_thkn = depth(2:end) - depth(1:end-1);
depth_mid = depth(1:end-1) + depth_thkn/2;
aus8_coor.depth_thkn = depth_thkn;
aus8_coor.depth_mid = depth_mid;
save([data_path 'SACS_data/aus8_coor'], 'aus8_coor')

aus8_u_g.mean = NaN(length(lat_u), length(lon_u), length(depth_mid));
aus8_v_g.mean = NaN(length(lat_v), length(lon_v), length(depth_mid));
for t = 1 : 12
    aus8_u_g.(Months{t}) = NaN(length(lat_u),length(lon_u),...
        length(depth_mid));
    aus8_v_g.(Months{t}) = NaN(length(lat_v),length(lon_v),...
        length(depth_mid));
end

% get u at vertical halfway points
for m = 1 : length(lat_u)
    for n = 1 : length(lon_u)
        aus8_u_g.mean(m,n,:) = interp1(...
            depth,squeeze(aus8_u_g_dum.mean(m,n,:)),depth_mid);
        for t = 1 : 12
            aus8_u_g.(Months{t})(m,n,:) = interp1(...
                depth,squeeze(aus8_u_g_dum.(Months{t})(m,n,:)),...
                depth_mid);
        end
    end
    fprintf('lat_u = %7.3f \n', lat_u(m))
end

% get v at vertical halfway points
for m = 1 : length(lat_v)
    for n = 1 : length(lon_v)
        aus8_v_g.mean(m,n,:) = interp1(...
            depth,squeeze(aus8_v_g_dum.mean(m,n,:)),depth_mid);
        for t = 1 : 12
            aus8_v_g.(Months{t})(m,n,:) = interp1(...
                depth,squeeze(aus8_v_g_dum.(Months{t})(m,n,:)),...
                depth_mid);
        end
    end
    fprintf('lat_v = %7.3f \n', lat_v(m))
end


%% save
save([data_path 'SACS_data/aus8_u_0'], 'aus8_u_0')
save([data_path 'SACS_data/aus8_v_0'], 'aus8_v_0')
save([data_path 'SACS_data/aus8_u_0_botm'], 'aus8_u_0_botm')
save([data_path 'SACS_data/aus8_v_0_botm'], 'aus8_v_0_botm')
save([data_path 'SACS_data/aus8_u_g'], 'aus8_u_g')
save([data_path 'SACS_data/aus8_v_g'], 'aus8_v_g')

disp(['aus8_u_0 aus8_v_0 aus8_u_0_botm aus8_v_0_botm ' ...
    'aus8_u_g aus8_v_g DONE'])

