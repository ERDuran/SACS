%% calculate adjustment velocities
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_u_g'])
load([data_path 'SACS_data/aus8_v_g'])
load([data_path 'SACS_data/aus8_U_d'])
load([data_path 'SACS_data/aus8_V_d'])
load([data_path 'SACS_data/aus8_U_ek'])
load([data_path 'SACS_data/aus8_V_ek'])
load([data_path 'SACS_data/aus8_F'])
load([data_path 'SACS_data/aus8_Lap_phi'])
load([data_path 'SACS_data/aus8_U_g_itgr'])
load([data_path 'SACS_data/aus8_V_g_itgr'])

a = aus8_coor.a;
pi180 = aus8_coor.pi180;
depth = aus8_coor.depth;
depth_mid = aus8_coor.depth_mid;
depth_thkn = aus8_coor.depth_thkn;
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;
u_bottom = -aus8_coor.u_bottom;
v_bottom = -aus8_coor.v_bottom;
Months = aus8_coor.Months;
Months{13} = 'mean';


%% 
for t = 1 : length(Months)
    
    u_d_now = aus8_U_d.(Months{t}) ./ u_bottom;
    u_d_repmat = repmat(u_d_now,1,1,length(depth_mid));
    aus8_u_g_prime.(Months{t}) = aus8_u_g.(Months{t}) - u_d_repmat;
    
    v_d_now = aus8_V_d.(Months{t}) ./ v_bottom;
    v_d_repmat = repmat(v_d_now,1,1,length(depth_mid));
    aus8_v_g_prime.(Months{t}) = aus8_v_g.(Months{t}) - v_d_repmat;
end

% calculate vertical integration
lat_u_length = length(lat_u);
lon_u_length = length(lon_u);
U_g_prime = NaN(lat_u_length, lon_u_length);
lat_v_length = length(lat_v);
lon_v_length = length(lon_v);
V_g_prime = NaN(lat_v_length, lon_v_length);

for t = 1 : length(Months)
    % integrate U along z. U_z is in m^2/s
    for ii = 1 : lat_u_length
        for jj = 1 : lon_u_length
            u_mid_now = squeeze(aus8_u_g_prime.(Months{t})(ii,jj,:));
            u_mid_now_times_depth_thicknesses = ...
                u_mid_now .* depth_thkn;
            aus8_U_g_prime.(Months{t})(ii,jj) = ...
                nansum(u_mid_now_times_depth_thicknesses);
        end
    end
    
    % integrate V along z. V_z is in m^2/s
    for ii = 1 : lat_v_length
        for jj = 1 : lon_v_length
            v_mid_now = squeeze(aus8_v_g_prime.(Months{t})(ii,jj,:));
            v_mid_now_times_depth_thicknesses = ...
                v_mid_now .* depth_thkn;
            aus8_V_g_prime.(Months{t})(ii,jj) = ...
                nansum(v_mid_now_times_depth_thicknesses);
        end
    end
    
    
    % Re-add the Ekman drift
    aus8_U_prime.(Months{t}) = ...
        aus8_U_g_prime.(Months{t}) + aus8_U_ek.(Months{t});
    aus8_V_prime.(Months{t}) = ...
        aus8_V_g_prime.(Months{t}) + aus8_V_ek.(Months{t});
    
    aus8_U_prime.(Months{t})(aus8_coor.U_mask) = 0;
    aus8_V_prime.(Months{t})(aus8_coor.V_mask) = 0;
    % calculate divergence
    du = (aus8_U_prime.(Months{t})(:,2:end) - ...
        aus8_U_prime.(Months{t})(:,1:end-1));
    dx = NaN(size(du));
    for ii = 1 : length(lat_u)
        dx_now = a * cos(lat_u(ii) * pi180) .* ...
            (lon_u(2:end) - lon_u(1:end-1)) * pi180;
        dx(ii,:) = dx_now;
    end
    lat_v_repmat = repmat(lat_v,1,length(lon_v));
    dv = aus8_V_prime.(Months{t})(1:end-1,:).* ...
        cos(lat_v_repmat(1:end-1,:) * pi180) - ...
        aus8_V_prime.(Months{t})(2:end,:).* ...
        cos(lat_v_repmat(2:end,:) * pi180);
    dy = NaN(size(dv));
    for jj = 1 : length(lon_v)
        dy_now = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
        dy(:,jj) = dy_now;
    end
    lat_u_repmat = repmat(lat_u,1,length(lon_v));
    aus8_F_prime.(Months{t}) = ...
        du./dx + 1./cos(lat_u_repmat * pi180).*dv./dy;
%     aus8_F_prime.(Months{t})(aus8_coor.F_mask) = NaN;
    
    aus8_U_prime.(Months{t})(aus8_coor.U_mask) = NaN;
    aus8_V_prime.(Months{t})(aus8_coor.V_mask) = NaN;
    
    disp([Months{t} ' OK!'])
end


%%
fig_n = 1;
rowcols = [2 1];
rowcols_size = [16 10]; % cm
margs = [1 1 1 1]; % cm
gaps = [1 1]; % cm
cmaps_levels = 12;

for sp = 1 : rowcols(1)*rowcols(2)
    axis_setup{sp} = ...
        [lon_v(1) lon_v(end) ...
        lat_u(end) lat_u(1)];
    x{sp} = lon_v;
    y{sp} = lat_u;
    cmaps{sp} = flipud(othercolor('RdBu8', cmaps_levels));
end

data = {...
    aus8_F.mean-aus8_Lap_phi.mean, ...
    aus8_F_prime.mean};
data{1}(data{1}==0) = NaN;
data{2}(data{2}==0) = NaN;

minmax = {...
    [-10^-5 10^-5], ...
    [-10^-5 10^-5]};
titles = {...
    ['aus8 mean $\nabla\cdot\bf{V}-\nabla^{2}\phi$']
    ['aus8 mean $\nabla\cdot\bf{V''}$']};
font_size = 9;
fig_color = [0.7 0.7 0.7];

fig = pcolor_maker(...
    fig_n, rowcols, rowcols_size, margs, gaps, ...
    x, y, data, axis_setup, minmax, cmaps, ...
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


%%
save([data_path 'SACS_data/aus8_u_g_prime'], 'aus8_u_g_prime')
save([data_path 'SACS_data/aus8_v_g_prime'], 'aus8_v_g_prime')
save([data_path 'SACS_data/aus8_U_g_prime_itgr'], 'aus8_U_g_prime')
save([data_path 'SACS_data/aus8_V_g_prime_itgr'], 'aus8_V_g_prime')
save([data_path 'SACS_data/aus8_U_prime'], 'aus8_U_prime')
save([data_path 'SACS_data/aus8_V_prime'], 'aus8_V_prime')
disp(['aus8_u_g_prime aus8_v_g_prime aus8_U_g_prime aus8_V_g_prime ' ...
    'aus8_U_prime aus8_V_prime DONE'])

