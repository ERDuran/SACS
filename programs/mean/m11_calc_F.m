%% calculate horizontal divergence
% clc
path(pathdef)

% set up main directory
cd /home/z5100174/Desktop/MATLAB

% Add path to the data to be loaded
addpath(genpath('functions'))

% Add path to the data to be loaded
addpath cars_out

% Add path to the gsw TEOS-10 library (including subfolders)
addpath(genpath('teos10_library'))

clear 
% vn is a function that can store the name of a variable
% hence vn(x) returns 'x' as a string.
vn = @(x) inputname(1);

load aus8_ZD_method

a = aus8_ZD_method.a.value;
cst_lat = aus8_ZD_method.f.cst_lat;
pi180 = aus8_ZD_method.a.pi180;


%% calculate horizontal divergence 
U = aus8_ZD_method.U;
lat_u = aus8_ZD_method.lat_u;
lon_u = aus8_ZD_method.lon_u;

V = aus8_ZD_method.V;
lat_v = aus8_ZD_method.lat_v;
lon_v = aus8_ZD_method.lon_v;

U(isnan(U)) = 0;
V(isnan(V)) = 0;

du = (U(:,2:end) - U(:,1:end-1));

dx = NaN(size(du));
for ii = 1 : length(lat_u)
    dx_now = a * cos(lat_u(ii) * pi180) .* ...
        (lon_u(2:end) - lon_u(1:end-1)) * pi180;
    
    dx(ii,:) = dx_now;
end

lat_repmat = repmat(lat_v,1,length(lon_v));
dv = V(1:end-1,:).*cos(lat_repmat(1:end-1,:) * pi180) - ...
    V(2:end,:).*cos(lat_repmat(2:end,:) * pi180);

dy = NaN(size(dv));
for jj = 1 : length(lon_v)
    dy_now = a * (lat_v(1:end-1) - lat_v(2:end)) * pi180;
    
    dy(:,jj) = dy_now;
end

lat_u_repmat = repmat(lat_u,1,length(lon_v));
F = du./dx + 1./cos(lat_u_repmat * pi180).*dv./dy;


%%
aus8_ZD_method.lat_F = lat_u;
aus8_ZD_method.lon_F = lon_v;
aus8_ZD_method.F = F;

save('cars_out/aus8_ZD_method', 'aus8_ZD_method')
disp('aus8_ZD_method saved in cars_out/aus8_ZD_method !')

