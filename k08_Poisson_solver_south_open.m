%% Solve Poisson equation using iteration method
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/KDau_U'])
load([data_path 'SACS_data/KDau_V'])
load([data_path 'SACS_data/KDau_F'])


%%
a = aus8_coor.a;
pi180 = aus8_coor.pi180;
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;
U_mask_KDau = aus8_coor.U_mask_KDau;
V_mask_KDau = aus8_coor.V_mask_KDau;
F_mask_KDau = aus8_coor.F_mask_KDau;
Months = aus8_coor.Months;


%%
% apply south boundary condition
lat_phi = [lat_u; lat_u(end)-1/8];
lon_phi = lon_v;

% dx for U position
dx_U = NaN(length(lat_u), length(lon_u)-1);
for ii = 1 : length(lat_u)
    dx_now = a * cos(lat_u(ii) * pi180) .* ...
        (lon_u(2:end) - lon_u(1:end-1)) * pi180;
    
    dx_U(ii,:) = dx_now;
end


% dy for V position
dy_V = NaN(length(lat_v)-1, length(lon_v));
for jj = 1 : length(lon_v)
    dy_now = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
    
    dy_V(:,jj) = dy_now;
end

% include open boundary phi
% dx for phi position
dx_phi = NaN(length(lat_phi)-1, length(lon_phi)-1);
for ii = 1 : length(lat_phi)-1
    dx_now = a * cos(lat_phi(ii) * pi180) .* ...
        (lon_phi(2:end) - lon_phi(1:end-1)) * pi180;
    
    dx_phi(ii,:) = dx_now;
end


% dy for phi position
dy_phi = NaN(length(lat_phi)-1, length(lon_phi));
for jj = 1 : length(lon_phi)
    dy_now = a * (lat_phi(1:end-1) - lat_phi(2:end)) * pi180;
    
    dy_phi(:,jj) = dy_now;
end

lat_V_repmat = repmat(lat_v,1,length(lon_phi));
lat_U_repmat = repmat(lat_u,1,length(lon_phi));


%%
% initial conditions
K = 1;
% time step
delta = max([dy_phi(:); dx_phi(:)]);
Dt = 0.5 * delta^2 / K *10^-0.75;

% Initial Conditions
KDau_phi.mean = KDau_F.mean;

% Boundary Condition
KDau_phi.mean(end+1,:) = 0;

for t = 1 : 12
    KDau_phi.(Months{t}) = KDau_F.(Months{t});
    KDau_phi.(Months{t})(end+1,:) = 0;
end

% number of iterations
n_5000s = 0 : 5000 : 5000000;
tol = 0.0001;
n_iter = 0;

Months{13} = 'mean';

for t = 1 : length(Months)
    rel_error = 1;
    tic
    while rel_error > tol
        n_iter = n_iter + 1;
        phi_n = KDau_phi.(Months{t});
        
        % Calculate grad(phi) = U_d + V_d
        KDau_U_d.(Months{t}) = zeros(size(KDau_U.(Months{t})));
        KDau_V_d.(Months{t}) = zeros(size(KDau_V.(Months{t})));
        
        phi_n_x = phi_n(1:end-1, 2:end) - phi_n(1:end-1, 1:end-1);
        phi_n_y = phi_n(1:end-1, :) - phi_n(2:end, :);
        
        KDau_U_d.(Months{t})(:,2:end-1) = (phi_n_x) ./ dx_phi;
        KDau_V_d.(Months{t})(2:end,:) = (phi_n_y) ./ dy_phi;
        
        KDau_U_d.(Months{t})(U_mask_KDau) = 0;
        KDau_V_d.(Months{t})(V_mask_KDau) = 0;
        
        % Calculate Laplacian(phi) = div(U_d,V_d)
        dU = KDau_U_d.(Months{t})(:,2:end) - ...
            KDau_U_d.(Months{t})(:,1:end-1);
        
        dV = KDau_V_d.(Months{t})(1:end-1,:).* ...
            cos(lat_V_repmat(1:end-1,:) * pi180) - ...
            KDau_V_d.(Months{t})(2:end,:).* ...
            cos(lat_V_repmat(2:end,:) * pi180);
        
        KDau_Lap_phi.(Months{t}) = ...
            dU./dx_U + 1./cos(lat_U_repmat * pi180).*dV./dy_V;
        
        K_Lap_phi_minus_F = ...
            K * KDau_Lap_phi.(Months{t}) - KDau_F.(Months{t});
        
        KDau_phi.(Months{t}) = ...
            phi_n(1:end-1,:) + Dt * (K_Lap_phi_minus_F);
        KDau_phi.(Months{t})(end+1,:) = 0;
        
        % calculate relative error
        rel_error = ...
            max(max(abs(KDau_phi.(Months{t}) - phi_n))) / ...
            max(max(abs(KDau_phi.(Months{t}))));
        
        if find(ismember(n_5000s,n_iter))
            fprintf('Iterations number = %.f Tolerance = %.9f \n', ...
                n_iter, rel_error)
        end
    end
    toc
end


%%
KDau_F_mean_plot = KDau_F.mean;
KDau_F_mean_plot(F_mask_KDau) = NaN;
KDau_Lap_phi_mean_plot = KDau_Lap_phi.mean;
KDau_Lap_phi_mean_plot(F_mask_KDau) = NaN;
K_Lap_phi_minus_F_plot = K_Lap_phi_minus_F;
K_Lap_phi_minus_F_plot(F_mask_KDau) = NaN;
KDau_U_d_mean_plot = KDau_U_d.mean;
KDau_U_d_mean_plot(U_mask_KDau) = NaN;
KDau_V_d_mean_plot = KDau_V_d.mean;
KDau_V_d_mean_plot(V_mask_KDau) = NaN;

fig_n = 1;
rowcols = [2 3];
rowcols_size = [8 5]; % cm
margs = [1 1 1 1]; % cm
gaps = [1 1]; % cm
cmaps_levels = 12;
x_chc = {lon_u, lon_v};
x_ind = [2 2 2 1 2 1];
y_chc = {lat_u, lat_v};
y_ind = [1 1 1 1 2 2];

for sp = 1 : rowcols(1)*rowcols(2)
    axis_setup{sp} = ...
        [aus8_coor.lon(1) aus8_coor.lon(end) ...
        aus8_coor.lat(end) aus8_coor.lat(1)];
    x{sp} = x_chc{x_ind(sp)};
    y{sp} = y_chc{y_ind(sp)};
    cmaps{sp} = flipud(othercolor('RdBu8', cmaps_levels));
end
    
data = {...
    KDau_F_mean_plot, ...
    KDau_Lap_phi_mean_plot, ...
    K_Lap_phi_minus_F_plot, ...
    KDau_U_d_mean_plot, ...
    KDau_V_d_mean_plot, ...
    NaN(length(lat_v), length(lon_u))};
minmax = {...
    [-10^-5 10^-5], ...
    [-10^-5 10^-5], ...
    [-10^-5 10^-5], ...
    [-10 10], ...
    [-10 10], ...
    [0 0]};
titles = {...
    ['KDau mean $\nabla\cdot\bf{V}$'], ...
    ['KDau mean $\nabla^{2}\phi$'], ...
    ['KDau mean $\nabla^{2}\phi-\nabla\cdot\bf{V}$'], ...
    ['KDau mean $U_{d}$ $(m^{2}/s)$'], ...
    ['KDau mean $V_{d}$ $(m^{2}/s)$'], ...
    [' ']};
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


%% save stuff
save([data_path 'SACS_data/KDau_Lap_phi'], 'KDau_Lap_phi')
save([data_path 'SACS_data/KDau_U_d'], 'KDau_U_d')
save([data_path 'SACS_data/KDau_V_d'], 'KDau_V_d')
disp('KDau_F KDau_Lap_phi KDau_U_d KDau_V_d DONE')

