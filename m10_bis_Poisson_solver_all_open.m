%% Solve Poisson equation using iteration method
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_U'])
load([data_path 'SACS_data/aus8_V'])
load([data_path 'SACS_data/aus8_F'])


%%
a = aus8_coor.a;
pi180 = aus8_coor.pi180;
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;
U_mask = aus8_coor.U_mask;
V_mask = aus8_coor.V_mask;
F_mask = aus8_coor.F_mask;
Months = aus8_coor.Months;


%%
% apply open ocean boundary everywhere
lat_phi = [lat_u(1)+1/8; lat_u; lat_u(end)-1/8];
lon_phi = [lon_v(1)-1/8, lon_v, lon_v(end)+1/8];

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
dx_phi = NaN(length(lat_phi)-2, length(lon_phi)-1);
for ii = 2 : length(lat_phi)-1
    dx_now = a * cos(lat_phi(ii) * pi180) .* ...
        (lon_phi(2:end) - lon_phi(1:end-1)) * pi180;
    
    dx_phi(ii-1,:) = dx_now;
end


% dy for phi position
dy_phi = NaN(length(lat_phi)-1, length(lon_phi)-2);
for jj = 2 : length(lon_phi)-1
    dy_now = a * (lat_phi(1:end-1) - lat_phi(2:end)) * pi180;
    
    dy_phi(:,jj-1) = dy_now;
end

lat_V_repmat = repmat(lat_v,1,length(lon_phi)-2);
lat_U_repmat = repmat(lat_u,1,length(lon_phi)-2);


%%
% initial conditions
K = 1;
% time step
delta = max([dy_phi(:); dx_phi(:)]);
Dt = 0.5 * delta^2 / K *10^-0.75;

% zero grid (boundary condition)
aus8_phi.mean = zeros(length(lat_phi), length(lon_phi));

for t = 1 : 12
    aus8_phi.(Months{t}) = aus8_F.(Months{t});
    aus8_phi.(Months{t})(end+1,:) = 0;
end

% number of iterations
n_5000s = 0 : 5000 : 5000000;
tol = 1 * 10^-6;
n_iter = 0;

Months{13} = 'mean';

for t = 13%1 : length(Months)
    rel_error = 1;
    tic
    while rel_error > tol
        n_iter = n_iter + 1;
        phi_n = aus8_phi.(Months{t});
        
        % Calculate grad(phi) = U_d + V_d
        aus8_U_d.(Months{t}) = zeros(size(aus8_U.(Months{t})));
        aus8_V_d.(Months{t}) = zeros(size(aus8_V.(Months{t})));
        
        phi_n_x = phi_n(2:end-1, 2:end) - phi_n(2:end-1, 1:end-1);
        phi_n_y = phi_n(1:end-1, 2:end-1) - phi_n(2:end, 2:end-1);
        
        aus8_U_d.(Months{t}) = (phi_n_x) ./ dx_phi;
        aus8_V_d.(Months{t}) = (phi_n_y) ./ dy_phi;
        
        aus8_U_d.(Months{t})(U_mask) = 0;
        aus8_V_d.(Months{t})(V_mask) = 0;
        
        % Calculate Laplacian(phi) = div(U_d,V_d)
        dU = aus8_U_d.(Months{t})(:,2:end) - ...
            aus8_U_d.(Months{t})(:,1:end-1);
        
        dV = aus8_V_d.(Months{t})(1:end-1,:).* ...
            cos(lat_V_repmat(1:end-1,:) * pi180) - ...
            aus8_V_d.(Months{t})(2:end,:).* ...
            cos(lat_V_repmat(2:end,:) * pi180);
        
        aus8_Lap_phi.(Months{t}) = ...
            dU./dx_U + 1./cos(lat_U_repmat * pi180).*dV./dy_V;
        
        K_Lap_phi_minus_F = ...
            K * aus8_Lap_phi.(Months{t}) - aus8_F.(Months{t});
        
        aus8_phi.(Months{t}) = zeros(length(lat_phi), length(lon_phi));
        aus8_phi.(Months{t})(2:end-1,2:end-1) = ...
            phi_n(2:end-1,2:end-1) + Dt * (K_Lap_phi_minus_F);
        
        rel_error = max(max(abs(K_Lap_phi_minus_F)));
        
        if find(ismember(n_5000s,n_iter))
            fprintf('Iterations number = %.f Tolerance = %.9f \n', ...
                n_iter, rel_error)
        end
    end
    toc
end


%%
aus8_F_mean_plot = aus8_F.mean;
aus8_F_mean_plot(F_mask) = NaN;
aus8_Lap_phi_mean_plot = aus8_Lap_phi.mean;
aus8_Lap_phi_mean_plot(F_mask) = NaN;
K_Lap_phi_minus_F_plot = K_Lap_phi_minus_F;
K_Lap_phi_minus_F_plot(F_mask) = NaN;
aus8_U_d_mean_plot = aus8_U_d.mean;
aus8_U_d_mean_plot(U_mask) = NaN;
aus8_V_d_mean_plot = aus8_V_d.mean;
aus8_V_d_mean_plot(V_mask) = NaN;

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
    aus8_F_mean_plot, ...
    aus8_Lap_phi_mean_plot, ...
    K_Lap_phi_minus_F_plot, ...
    aus8_U_d_mean_plot, ...
    aus8_V_d_mean_plot, ...
    NaN(length(lat_v), length(lon_u))};
minmax = {...
    [-10^-5 10^-5], ...
    [-10^-5 10^-5], ...
    [-10^-7 10^-7], ...
    [-10 10], ...
    [-10 10], ...
    [0 0]};
titles = {...
    ['aus8 mean $\nabla\cdot\bf{V}$'], ...
    ['aus8 mean $\nabla^{2}\phi$'], ...
    ['aus8 mean $\nabla^{2}\phi-\nabla\cdot\bf{V}$'], ...
    ['aus8 mean $U_{d}$ $(m^{2}/s)$'], ...
    ['aus8 mean $V_{d}$ $(m^{2}/s)$'], ...
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
save([data_path 'SACS_data/aus8_Lap_phi'], 'aus8_Lap_phi')
save([data_path 'SACS_data/aus8_U_d'], 'aus8_U_d')
save([data_path 'SACS_data/aus8_V_d'], 'aus8_V_d')
disp('aus8_F aus8_Lap_phi aus8_U_d aus8_V_d DONE')

