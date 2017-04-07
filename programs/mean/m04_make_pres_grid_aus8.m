%% project data in aus8 coast and a8cst TS into hr pressure grid
clc
path(pathdef)

% set up main directory
cd /home/z5100174/Desktop/MATLAB

% Add path to the data to be loaded
addpath(genpath('functions'))

% Add path to the data to be loaded
addpath cars_out

% Add path to the gsw TEOS-10 library (including subfolders)
addpath(genpath('teos10_library'))


%% Load a8c
clear 
load aus8, disp('aus8 loaded !')

data_names = fieldnames(aus8);


%% project ac8 into hr pres grid
for a = 1 : length(data_names)
        lat = aus8.(data_names{a}).var.lat;
        lon = aus8.(data_names{a}).var.lon;
        depth = aus8.(data_names{a}).var.depth;
        mean_raw = aus8.(data_names{a}).var.mean;
        
        pres = (0 : 10 : 2000)';
        
        mean_on_pres = NaN(length(lat), length(lon), length(pres));
        
        for ii = 1 : length(lat)
            depth_on_pres = -gsw_z_from_p(pres, lat(ii));
            
            for jj = 1 : length(lon)
                 mean_on_pres(ii,jj,:) = ...
                     interp1(depth, squeeze(mean_raw(ii,jj,:)), depth_on_pres);
                
            end
        end
        
        aus8.(data_names{a}).var_on_pres.lat = lat;
        aus8.(data_names{a}).var_on_pres.lon = lon;
        aus8.(data_names{a}).var_on_pres.pres = pres;
        aus8.(data_names{a}).var_on_pres.mean = mean_on_pres;
        fprintf('%s projected onto pressure grid !\n', data_names{a})
        
end


%% save updated a8c
save('cars_out/aus8', 'aus8')
disp('CARS aus8 saved into aus8 !')

