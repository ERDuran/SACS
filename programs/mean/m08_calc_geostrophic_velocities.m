%% calculate geostrophic velocity assuming f is constant
clearvars('-except', '*_path')

load aus8_streamfunction
lat_v = aus8_streamfunction.lat;
lon_u = aus8_streamfunction.lon;
pres = aus8_streamfunction.pres;
p_ref_list = aus8_streamfunction.p_ref_list;
aus8_geostrophy.pres = pres;
a = 6371000; % Earth's radius in meters
pi180 = pi/180;


%% set up the parameters
aus8_geostrophy.a.name = 'Earth radius';
aus8_geostrophy.a.symbol = 'a';
aus8_geostrophy.a.units = 'm';
aus8_geostrophy.a.pi180 = pi180;
aus8_geostrophy.a.value = a;

% mid lat grid
lat_u = (lat_v(1:end-1) + lat_v(2:end))/2;

% mid lon grid
lon_v = (lon_u(1:end-1) + lon_u(2:end))/2;

aus8_geostrophy.u.name = 'Zonal Geostrophic Velocity';
aus8_geostrophy.u.symbol = 'u_g';
aus8_geostrophy.u.units = 'cm/s';
aus8_geostrophy.u.lat_u = lat_u;
aus8_geostrophy.u.lon_u = lon_u;

aus8_geostrophy.v.name = 'Meridional Geostrophic Velocity';
aus8_geostrophy.v.symbol = 'v_g';
aus8_geostrophy.v.units = 'cm/s';
aus8_geostrophy.v.lat_v = lat_v;
aus8_geostrophy.v.lon_v = lon_v;

dynh_raw = aus8_streamfunction.dynh.p_ref_0;
dynh = dynh_raw;

% smooth out dynh
for ii = 2 : length(lat_v)-1
    for jj = 2 : length(lon_u)-1
        finite_dynh = find(isfinite(dynh_raw(ii,jj,:)));
        
        if ~isempty(finite_dynh)
            %
            dynh_raw_N = squeeze(dynh_raw(ii-1,jj,finite_dynh));
            dynh_raw_S = squeeze(dynh_raw(ii+1,jj,finite_dynh));
            dynh_raw_W = squeeze(dynh_raw(ii,jj-1,finite_dynh));
            dynh_raw_E = squeeze(dynh_raw(ii,jj+1,finite_dynh));
            dynh_raw_C = squeeze(dynh_raw(ii,jj,finite_dynh));
            
            %
            % %%% 1 way: if any of the four points around are nan then
            % the centre point is the original value ie. only averaging
            % over the 5 points
            % PRETTY GOOD WAY
            %
            % mean_dynh = mean([...
            %     dynh_raw_N, ...
            %     dynh_raw_S, ...
            %     dynh_raw_W, ...
            %     dynh_raw_E, ...
            %     dynh_raw_C],2);
            %
            % mean_dynh_NaN = isnan(mean_dynh);
            % mean_dynh(mean_dynh_NaN) = dynh_raw_C(mean_dynh_NaN);
            % dynh(ii,jj,finite_dynh) = mean_dynh;
            % %%%
            %
            %
            % %%% 2 way: average along the whole finite water column of
            % % the centre ie. allowed to average over less than 5
            % % points BAD WAY
            % mean_dynh = nanmean([...
            %     dynh_raw_N, ...
            %     dynh_raw_S, ...
            %     dynh_raw_W, ...
            %     dynh_raw_E, ...
            %     dynh_raw_C],2);
            %
            % %
            % dynh(ii,jj,finite_dynh) = mean_dynh;
            % %%%
            
            
            %%% 3 way: if we are above the shelf (ie. depth shallower
            % than 2000 db), then do not average at all and take
            % original centre cell. the average offshore is always over
            % 5 points
            % SLIGHTLY BETTER THAN 1
            if finite_dynh(end) == length(pres)
                %
                mean_dynh = mean([...
                    dynh_raw_N, ...
                    dynh_raw_S, ...
                    dynh_raw_W, ...
                    dynh_raw_E, ...
                    dynh_raw_C],2);
                
                %
                mean_dynh_NaN = isnan(mean_dynh);
                mean_dynh(mean_dynh_NaN) = dynh_raw_C(mean_dynh_NaN);
                dynh(ii,jj,finite_dynh) = mean_dynh;
                
            else
                %
                dynh(ii,jj,finite_dynh) = dynh_raw_C;
                
            end
            %%%
            
            %
            % %%% 4 way: if we are above the shelf (ie. depth shallower
            % % than 2000 db), then do not average at all and take
            % % original centre cell. the average offshore can be over
            % % less than 5 points
            % % PRETTY BAD WAY
            % if finite_dynh(end) == length(pres)
            %     %
            %     mean_dynh = nanmean([...
            %         dynh_raw_N, ...
            %         dynh_raw_S, ...
            %         dynh_raw_W, ...
            %         dynh_raw_E, ...
            %         dynh_raw_C],2);
            %
            %     %
            %     dynh(ii,jj,finite_dynh) = mean_dynh;
            %
            % else
            %     %
            %     dynh(ii,jj,finite_dynh) = dynh_raw_C;
            %
            % end
            % %%%
        end
    end
end


%% Calculate geostrophic velocities
% u
disp('Calculating geostrophic velocities ...')
u_na = NaN(length(lat_v)-1,length(lon_u),length(pres));
f_u = -gsw_f(lat_u);
f_u_repmat = repmat(f_u, 1, length(pres));
for jj = 1 : length(lon_u)
    % dynamic height at that longitude
    dynh_now = squeeze(dynh(:,jj,:));
    
    % dynamic height difference
    dynh_diff_now = ...
        dynh_now(1:end-1,:) - dynh_now(2:end,:);
    
    % latitude distance
    dy_now = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
    dy_now_repmat = repmat(dy_now,1,length(pres));
    
    % u derivation
    u_now = dynh_diff_now ./ (dy_now_repmat .* f_u_repmat);
    
    % save the current u and turn into cm/s
    u_na(:,jj,:) = u_now;
    fprintf('lon = %7.3f \n', lon_u(jj))
    
end
% save u
aus8_geostrophy.u.u_0 = u_na;


% v
v_na = NaN(length(lat_v),length(lon_u)-1,length(pres));
for ii = 1 : length(lat_v)
    f_v = gsw_f(lat_v(ii));
    f_v_repmat = repmat(f_v, length(lon_v), length(pres));

    % dynamic height at that latitude
    dynh_now = squeeze(dynh(ii,:,:));
    
    % dynamic height difference
    dynh_diff_now = ...
        dynh_now(2:end,:) - dynh_now(1:end-1,:);
    
    % longitude distance
    dx_now = a * cos(lat_v(ii) * pi180) * ...
        (lon_u(2:end) - lon_u(1:end-1)) * pi180;
    dx_now_repmat = repmat(dx_now',1,length(pres));
    
    % v derivation
    v_now = dynh_diff_now ./ (dx_now_repmat .* f_v_repmat);
    
    % save the current v and convert into cm/s
    v_na(ii,:,:) = v_now;
    fprintf('lat = %7.3f \n', lat_v(ii))
    
end
% save v
aus8_geostrophy.v.v_0 = v_na;
fprintf('Geostrophic Velocity Calculated')


% Adjust the geostrophic velocities
u_barotropic = NaN(size(u_na));
for ii = 1 : length(lat_u)
    for jj = 1 : length(lon_u)
        last_finite = find(isfinite(u_na(ii,jj,:)),1,'last');
        
        if isempty(last_finite)
            continue
            
        else
            u_bottom = u_na(ii,jj,last_finite);
            u_barotropic(ii,jj,:) = repmat(u_bottom,length(pres),1);
            
        end
    end
end
u_HH = u_na - u_barotropic;

v_barotropic = NaN(size(v_na));
for ii = 1 : length(lat_v)
    for jj = 1 : length(lon_v)
        last_finite = find(isfinite(v_na(ii,jj,:)),1,'last');
        
        if isempty(last_finite)
            continue
            
        else
            v_bottom = v_na(ii,jj,last_finite);
            v_barotropic(ii,jj,:) = repmat(v_bottom,length(pres),1);
            
        end
    end
end
v_HH = v_na - v_barotropic;

aus8_geostrophy.u.u_0_barotropic_bottom_2000 = ...
    u_barotropic(:,:,1);
aus8_geostrophy.u.u_0_HH_2000 = u_HH;

aus8_geostrophy.v.v_0_barotropic_bottom_2000 = ...
    v_barotropic(:,:,1);
aus8_geostrophy.v.v_0_HH_2000 = v_HH;


% save
save([cars_out_path 'aus8_geostrophy'], 'aus8_geostrophy')
disp(['aus8_geostrophy saved in ' ...
    cars_out_path 'aus8_geostrophy'])

