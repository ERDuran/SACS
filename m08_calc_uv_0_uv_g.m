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
depth_thkn = -(depth(2:end) - depth(1:end-1));
depth_mid = depth(1:end-1) - depth_thkn/2;
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


%% figure set-up
fig_n = 1;
rowcols = [3 3];
rowcols_size = [8 5]; % cm
margs = [1 1 1 1]; % cm
gaps = [1 1]; % cm

z = {0 -400 -1000 0 -400 -1000 0 -400 -1000};
for p = 1 : length(z)
    z_ind{p} = find_nearest(aus8_coor.depth_mid,z{p});
end
cmaps_levels = 12;

cmaps_chc = {'RdBu8', 'BuPu8'};
cmaps_ind = [1 1 1 1 1 1 2 2 2];
x_chc = {lon_u, lon_v, lon};
x_ind = [1 1 1 2 2 2 3 3 3];
y_chc = {lat_u, lat_v, lat};
y_ind = [1 1 1 2 2 2 3 3 3];
minmax_chc = {...
    [-0.1 0.1], ...
    [-0.1 0.1], ...
    [0 0.2]};

for sp = 1 : rowcols(1)*rowcols(2)
    axis_setup{sp} = ...
        [lon(1) lon(end) ...
        lat(end) lat(1)];
    cmaps{sp} = ...
        othercolor(cmaps_chc{cmaps_ind(sp)}, cmaps_levels);
    if cmaps_ind(sp) == 1
        cmaps{sp} = flipud(cmaps{sp});
    end
    x{sp} = x_chc{x_ind(sp)};
    y{sp} = y_chc{y_ind(sp)};
    minmax{sp} = minmax_chc{x_ind(sp)};
    nn{sp} = 12;
    s{sp} = 2;
end

% data sorting
for sp = 1 : rowcols(1)*rowcols(2)
    if cmaps_ind(sp) == 1
        if x_ind(sp) == 1
            u{sp} = NaN(size(aus8_u_g.mean(:,:,1)));
            v{sp} = NaN(size(aus8_u_g.mean(:,:,1)));
            data{sp} = aus8_u_g.mean(:,:,z_ind{sp});
            titles{sp} = ...
                ['aus8 mean HH $u_{g}$ $z=' ...
                num2str(aus8_coor.depth_mid(z_ind{sp})) '$ $(m/s)$'];

        elseif x_ind(sp) == 2
            u{sp} = NaN(size(aus8_v_g.mean(:,:,1)));
            v{sp} = NaN(size(aus8_v_g.mean(:,:,1)));
            data{sp} = aus8_v_g.mean(:,:,z_ind{sp});
            titles{sp} = ...
                ['aus8 mean HH $v_{g}$ $z=' ...
                num2str(aus8_coor.depth_mid(z_ind{sp})) '$ $(m/s)$'];
        end

    elseif cmaps_ind(sp) == 2
        u{sp} = interp2(lon_u,lat_u, ...
                aus8_u_g.mean(:,:,z_ind{sp}),lon,lat);
        v{sp} = interp2(lon_v,lat_v, ...
                aus8_v_g.mean(:,:,z_ind{sp}),lon,lat);
        data{sp} = sqrt(u{sp}.^2 + v{sp}.^2);
        titles{sp} = ...
            ['aus8 mean HH $\bf{v_{g}}$ $z=' ...
            num2str(aus8_coor.depth_mid(z_ind{sp})) '$ $(m/s)$'];
    end

end

font_size = 9;
fig_color = [0.7 0.7 0.7];

fig = quiver_maker(...
    fig_n, rowcols, rowcols_size, margs, gaps, ...
    x, y, data, u, v, nn, s, axis_setup, minmax, cmaps, ...
    titles, font_size, fig_color);
outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
export_fig(fig, ...
    [figures_path mfilename '/' scriptname(1:3) ...
    '_fig' num2str(fig_n) '_'], ...
    '-m3')
close


%% figure set-up
fig_n = 2;
rowcols = [3 1];
rowcols_size = [25 9.5]; % cm
margs = [1 1 1 1]; % cm
gaps = [1 1]; % cm

z = {0 -400 -1000 0 -400 -1000 0 -400 -1000};
for p = 1 : length(z)
    z_ind{p} = find_nearest(aus8_coor.depth_mid,z{p});
end
cmaps_levels = 12;

cmaps_chc = {'RdBu8', 'BuPu8'};
cmaps_ind = [2 2 2];
x_chc = {lon_u, lon_v, lon};
x_ind = [3 3 3];
y_chc = {lat_u, lat_v, lat};
y_ind = [3 3 3];
minmax_chc = {...
    [-0.1 0.1], ...
    [-0.1 0.1], ...
    [0 0.2]};

for sp = 1 : rowcols(1)*rowcols(2)
    axis_setup{sp} = ...
        [113 148 ...
        -47 -31.5];
    cmaps{sp} = ...
        othercolor(cmaps_chc{cmaps_ind(sp)}, cmaps_levels);
    if cmaps_ind(sp) == 1
        cmaps{sp} = flipud(cmaps{sp});
    end
    x{sp} = x_chc{x_ind(sp)};
    y{sp} = y_chc{y_ind(sp)};
    minmax{sp} = minmax_chc{x_ind(sp)};
end

% data sorting
for sp = 1 : rowcols(1)*rowcols(2)
    if cmaps_ind(sp) == 1
        if x_ind(sp) == 1
            u{sp} = NaN(size(aus8_u_g.mean(:,:,1)));
            v{sp} = NaN(size(aus8_u_g.mean(:,:,1)));
            data{sp} = aus8_u_g.mean(:,:,z_ind{sp});
            titles{sp} = ...
                ['aus8 mean HH $u_{g}$ $z=' ...
                num2str(aus8_coor.depth_mid(z_ind{sp})) '$ $(m/s)$'];

        elseif x_ind(sp) == 2
            u{sp} = NaN(size(aus8_v_g.mean(:,:,1)));
            v{sp} = NaN(size(aus8_v_g.mean(:,:,1)));
            data{sp} = aus8_v_g.mean(:,:,z_ind{sp});
            titles{sp} = ...
                ['aus8 mean HH $v_{g}$ $z=' ...
                num2str(aus8_coor.depth_mid(z_ind{sp})) '$ $(m/s)$'];
        end

    elseif cmaps_ind(sp) == 2
        u{sp} = interp2(lon_u,lat_u, ...
                aus8_u_g.mean(:,:,z_ind{sp}),lon,lat);
        v{sp} = interp2(lon_v,lat_v, ...
                aus8_v_g.mean(:,:,z_ind{sp}),lon,lat);
        data{sp} = sqrt(u{sp}.^2 + v{sp}.^2);
        titles{sp} = ...
            ['aus8 mean HH $\bf{v_{g}}$ $z=' ...
            num2str(aus8_coor.depth_mid(z_ind{sp})) '$ $(m/s)$'];
    end

end

nn = {3 4 4};
s = {12 7 6};

font_size = 9;
fig_color = [0.7 0.7 0.7];

fig = quiver_maker(...
    fig_n, rowcols, rowcols_size, margs, gaps, ...
    x, y, data, u, v, nn, s, axis_setup, minmax, cmaps, ...
    titles, font_size, fig_color);
outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
export_fig(fig, ...
    [figures_path mfilename '/' scriptname(1:3) ...
    '_fig' num2str(fig_n) '_'], ...
    '-m3')
close


%% save
save([data_path 'SACS_data/aus8_u_0'], 'aus8_u_0')
save([data_path 'SACS_data/aus8_v_0'], 'aus8_v_0')
save([data_path 'SACS_data/aus8_u_0_botm'], 'aus8_u_0_botm')
save([data_path 'SACS_data/aus8_v_0_botm'], 'aus8_v_0_botm')
save([data_path 'SACS_data/aus8_u_g'], 'aus8_u_g')
save([data_path 'SACS_data/aus8_v_g'], 'aus8_v_g')

disp(['aus8_u_0 aus8_v_0 aus8_u_0_botm aus8_v_0_botm ' ...
    'aus8_u_g aus8_v_g DONE'])

