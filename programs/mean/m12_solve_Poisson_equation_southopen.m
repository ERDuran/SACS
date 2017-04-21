%% Solve Poisson equation using iteration method
clearvars('-except', '*_path')

load aus8_ZD_method
a = aus8_ZD_method.a;
pi180 = aus8_ZD_method.pi180;
lat_U = aus8_ZD_method.lat_u;
lon_U = aus8_ZD_method.lon_u;
U = aus8_ZD_method.U;
lat_V = aus8_ZD_method.lat_v;
lon_V = aus8_ZD_method.lon_v;
V = aus8_ZD_method.V;
F = aus8_ZD_method.F;


%%
U_nan = isnan(U);
U(U_nan) = 0;
V_nan = isnan(V);
V(V_nan) = 0;

% apply south boundary condition
lat_phi = [lat_U; lat_U(end)-1/8];
lon_phi = lon_V;

%
F_nan = F == 0;

% dx for U position
dx_U = NaN(length(lat_U), length(lon_U)-1);
for ii = 1 : length(lat_U)
    dx_now = a * cos(lat_U(ii) * pi180) .* ...
        (lon_U(2:end) - lon_U(1:end-1)) * pi180;
    
    dx_U(ii,:) = dx_now;
end


% dy for V position
dy_V = NaN(length(lat_V)-1, length(lon_V));
for jj = 1 : length(lon_V)
    dy_now = a * (lat_V(1:end-1) - lat_V(2:end)) * pi180;
    
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

lat_V_repmat = repmat(lat_V,1,length(lon_phi));
lat_U_repmat = repmat(lat_U,1,length(lon_phi));


%%
% initial conditions
K = 1;
% time step
delta = max([dy_phi(:); dx_phi(:)]);
Dt = 0.5 * delta^2 / K *10^-0.75;

% Initial Conditions
phi_np1 = ones(size(F));
phi_np1(isnan(F)) = NaN;

% Boundary Condition
phi_np1(end+1,:) = 0;

% number of iterations
n = 100000;
n_10000s = 0 : 10000 : n;

rel_error = NaN(n,1);
max_diff_K_Lap_phi__F = NaN(n,1);
number_NaN = NaN(n,1);

tic
for nn = 1 : n
    phi_n = phi_np1;
    
    % Calculate grad(phi) = U_d + V_d
    U_d = zeros(size(U));
    V_d = zeros(size(V));
    
    phi_n_x = phi_n(1:end-1, 2:end) - phi_n(1:end-1, 1:end-1);
    phi_n_y = phi_n(1:end-1, :) - phi_n(2:end, :);
    
    U_d(:,2:end-1) = (phi_n_x) ./ dx_phi;
    V_d(2:end,:) = (phi_n_y) ./ dy_phi;
    
    U_d(U_nan) = 0;
    V_d(V_nan) = 0;
    
%     % manually add boundary condition on the domain edge (stencil
%     % does not scan over the edge)
%     U_d_dom(1,54) = 0; 
%     V_d_dom(16,end) = 0;
%     
%     U_d_dom(1,1) = 0;
%     U_d_dom(end,1) = 0;
%     U_d_dom(end,end) = 0;
%     
%     V_d_dom(1,1) = 0;
    
    % Calculate Laplacian(phi) = div(U_d,V_d)
    du = (U_d(:,2:end) - U_d(:,1:end-1));
    
    dv = V_d(1:end-1,:).*cos(lat_V_repmat(1:end-1,:) * pi180) - ...
        V_d(2:end,:).*cos(lat_V_repmat(2:end,:) * pi180);
    
    Lap_phi = du./dx_U + 1./cos(lat_U_repmat * pi180).*dv./dy_V;
       
    K_Lap_phi_minus_F = K * Lap_phi - F;
    
    phi_np1 = phi_n(1:end-1,:) + Dt * (K_Lap_phi_minus_F);
    phi_np1(end+1,:) = 0;
    
    number_NaN(nn) = length(find(isnan(phi_np1)));
    
    %
    max_diff_K_Lap_phi__F(nn) = max(max(K_Lap_phi_minus_F));
    
    % calculate relative error
    rel_error(nn) = ...
        max(max(abs(phi_np1 - phi_n))) / max(max(abs(phi_np1)));
        
    if find(ismember(n_10000s,nn))
        fprintf('iteration number = %7.0f \n', nn)
    end
end
toc

%
phi = phi_np1;

%
U_d(U_nan) = NaN;
V_d(V_nan) = NaN;

%
F(F_nan) = NaN;
Lap_phi(F_nan) = NaN;
K_Lap_phi_minus_F(F_nan) = NaN;


%%
close
figure
set(gcf,'units','normalized','outerposition',[0 0.05 0.95 0.95])

subplot(2,3,1)
m = 5;
[CMAP,LEV,WID,CLIM] = cmjoin({...
    flipud(othercolor('Blues7')),othercolor('Reds7')}, ...
    [-10^-m 0 10^-m], [10^(-m-1) 10^(-m-1)]);
colormap(CMAP)
pcolor(F), shading interp
colorbar, axis ij, title('F')
caxis(CLIM)

subplot(2,3,2)
colormap(CMAP)
pcolor(Lap_phi), shading interp
colorbar, axis ij, title('Lap(phi)')
caxis(CLIM)

subplot(2,3,3)
m = 6;
[CMAP,LEV,WID,CLIM] = cmjoin({...
    flipud(othercolor('Blues7')),othercolor('Reds7')}, ...
    [-10^-m 0 10^-m], [10^(-m-1) 10^(-m-1)]);
pcolor(K_Lap_phi_minus_F), shading interp
colorbar, axis ij, title('K * Lap(phi) - F')
caxis(CLIM)

subplot(2,3,4)
m = -1;
[CMAP,LEV,WID,CLIM] = cmjoin({...
    flipud(othercolor('Blues7')),othercolor('Reds7')}, ...
    [-10^-m 0 10^-m], [10^(-m-1) 10^(-m-1)]);
colormap(CMAP)
pcolor(U_d), shading interp
colorbar, axis ij, title('U_d')
caxis(CLIM)

subplot(2,3,5)
pcolor(V_d), shading interp
colorbar, axis ij, title('V_d')
caxis(CLIM)

subplot(2,3,6)
text(0.1,0.8,['Dt = ' num2str(Dt(end))])

text(0.1,0.6,['Number of iterations = ' num2str(nn)])

text(0.1,0.4,['Relative error = ' num2str(rel_error(end))])

% Save
outputls = ls(figures_path);
scriptname = mfilename;
if ~contains(outputls, scriptname)
    mkdir(figures_path, scriptname)
end
export_fig(fig1, [figures_path mfilename '/' scriptname(1:3) '_'], ...
    '-m4')
close


%% save stuff
aus8_ZD_method.lat_phi = lat_phi;
aus8_ZD_method.lon_phi = lon_phi;
aus8_ZD_method.phi = phi;
aus8_ZD_method.U_d = U_d;
aus8_ZD_method.V_d = V_d;
aus8_ZD_method.Lap_phi = Lap_phi;
aus8_ZD_method.n_iteration = nn;
aus8_ZD_method.rel_error = rel_error;
aus8_ZD_method.max_diff_K_Lap_phi__F = max_diff_K_Lap_phi__F;
aus8_ZD_method.Dt = Dt;
aus8_ZD_method.K = K;

save([cars_out_path 'aus8_ZD_method'], 'aus8_ZD_method')
disp(['aus8_ZD_method saved in ' ...
    cars_out_path 'aus8_ZD_method'])

