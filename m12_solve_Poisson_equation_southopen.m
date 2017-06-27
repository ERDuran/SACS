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
aus8_phi.mean = aus8_F.mean;

% Boundary Condition
aus8_phi.mean(end+1,:) = 0;

% number of iterations
n_5000s = 0 : 5000 : 5000000;
n_iter = 0;
tol = 0.0001;
rel_error = 1;

tic
while rel_error > tol
    n_iter = n_iter + 1;
    phi_n = aus8_phi.mean;
    
    % Calculate grad(phi) = U_d + V_d
    aus8_U_d.mean = zeros(size(aus8_U.mean));
    aus8_V_d.mean = zeros(size(aus8_V.mean));
    
    phi_n_x = phi_n(1:end-1, 2:end) - phi_n(1:end-1, 1:end-1);
    phi_n_y = phi_n(1:end-1, :) - phi_n(2:end, :);
    
    aus8_U_d.mean(:,2:end-1) = (phi_n_x) ./ dx_phi;
    aus8_V_d.mean(2:end,:) = (phi_n_y) ./ dy_phi;
    
    aus8_U_d.mean(U_mask) = 0;
    aus8_V_d.mean(V_mask) = 0;
    
    % Calculate Laplacian(phi) = div(U_d,V_d)
    dU = aus8_U_d.mean(:,2:end) - aus8_U_d.mean(:,1:end-1);
    
    dV = aus8_V_d.mean(1:end-1,:).* ...
        cos(lat_V_repmat(1:end-1,:) * pi180) - ...
        aus8_V_d.mean(2:end,:).* ...
        cos(lat_V_repmat(2:end,:) * pi180);
    
    aus8_Lap_phi.mean = dU./dx_U + 1./cos(lat_U_repmat * pi180).*dV./dy_V;
       
    K_Lap_phi_minus_F = K * aus8_Lap_phi.mean - aus8_F.mean;
    
    aus8_phi.mean = phi_n(1:end-1,:) + Dt * (K_Lap_phi_minus_F);
    aus8_phi.mean(end+1,:) = 0; 
        
    % calculate relative error
    rel_error = ...
        max(max(abs(aus8_phi.mean - phi_n))) / max(max(abs(aus8_phi.mean)));
        
    if find(ismember(n_5000s,n_iter))
        fprintf('Iterations number = %.f Tolerance = %.9f \n', ...
            n_iter, rel_error)
    end
end
toc


%
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

close
fig1 = figure;
set(gcf,'units','normalized','outerposition',[0 0.05 0.95 0.95])

subplot(2,3,1)
m = 5;
[CMAP,LEV,WID,CLIM] = cmjoin({...
    flipud(othercolor('Blues7')),othercolor('Reds7')}, ...
    [-10^-m 0 10^-m], [10^(-m-1) 10^(-m-1)]);
colormap(CMAP)
pcolor(aus8_F_mean_plot), shading interp
colorbar, axis ij, title('F')
caxis(CLIM)

subplot(2,3,2)
colormap(CMAP)
pcolor(aus8_Lap_phi_mean_plot), shading interp
colorbar, axis ij, title('Lap(phi)')
caxis(CLIM)

subplot(2,3,3)
[CMAP,LEV,WID,CLIM] = cmjoin({...
    flipud(othercolor('Blues7')),othercolor('Reds7')}, ...
    [-10^-m 0 10^-m], [10^(-m-1) 10^(-m-1)]);
pcolor(K_Lap_phi_minus_F_plot), shading interp
colorbar, axis ij, title('K * Lap(phi) - F')
caxis(CLIM)

subplot(2,3,4)
m = -1;
[CMAP,LEV,WID,CLIM] = cmjoin({...
    flipud(othercolor('Blues7')),othercolor('Reds7')}, ...
    [-10^-m 0 10^-m], [10^(-m-1) 10^(-m-1)]);
colormap(CMAP)
pcolor(aus8_U_d_mean_plot), shading interp
colorbar, axis ij, title('U_d')
caxis(CLIM)

subplot(2,3,5)
pcolor(aus8_V_d_mean_plot), shading interp
colorbar, axis ij, title('V_d')
caxis(CLIM)

subplot(2,3,6)
text(0.1,0.8,['Dt = ' num2str(Dt)])

text(0.1,0.6,['Number of iterations = ' num2str(n_iter)])

text(0.1,0.4,['Relative error = ' num2str(rel_error)])

% Save
% outputls = ls(figures_path);
% scriptname = mfilename;
% if ~contains(outputls, scriptname)
%     mkdir(figures_path, scriptname)
% end
% export_fig(fig1, [figures_path mfilename '/' scriptname(1:3) '_'], ...
%     '-m4')
% close


%% save stuff
aus8_ZD_method.lat_phi = lat_phi;
aus8_ZD_method.lon_phi = lon_phi;
aus8_ZD_method.phi = phi;
aus8_ZD_method.U_d = aus8_U_d;
aus8_ZD_method.V_d = aus8_V_d;
aus8_ZD_method.Lap_phi = aus8_Lap_phi;
aus8_ZD_method.n_iteration = n_iter;
aus8_ZD_method.rel_error = rel_error;
aus8_ZD_method.Dt = Dt;
aus8_ZD_method.K = K;

save([cars_out_path 'aus8_ZD_method'], 'aus8_ZD_method')
disp(['aus8_ZD_method saved in ' ...
    cars_out_path 'aus8_ZD_method'])

