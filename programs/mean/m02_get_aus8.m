%% Subset data off CARS aus8 at 3 locations around Australia's coast
clc
path(pathdef)

% set up main directory
cd /home/z5100174/Desktop/MATLAB

% Add path of useful functions (including subfolders)
addpath(genpath('functions'))

% Add path of the calculations output folder
addpath cars_out


%% Load aus8
clear
% vn is a function that can store the name of a variable
% hence vn(x) returns 'x' as a string.
vn = @(x) inputname(1);

load aus8_raw
disp('CARS aus8_raw loaded !')


%% projections, domains definition
% Lat subset
lat_subset = -30 : -1/8 : -50;
% Lon subset
lon_subset = 108 : 1/8 : 153;

% Depth subset
depth_limit = 2250;


%% subset from aus8 and into aus8_coast
aus8_names = fieldnames(aus8_raw);
var_names = aus8_raw.temp.var_names;

for a = 1 : length(aus8_names)
    lon_aus8 = aus8_raw.(aus8_names{a}).var.lon;
    lat_aus8 = aus8_raw.(aus8_names{a}).var.lat;
    depth_aus8 = aus8_raw.(aus8_names{a}).var.depth;
    
    % lon index
    lon_ind = ismember(lon_aus8, lon_subset);
    lon = lon_aus8(lon_ind);
    
    % lat index
    lat_ind = ismember(lat_aus8, lat_subset);
    lat = lat_aus8(lat_ind);
    
    % depth index
    depth_ind = depth_aus8 <= depth_limit;
    depth = depth_aus8(depth_ind);
    
    for c = 1 : length(var_names)
        if ismember(1, strcmp(var_names{c}, {'an_cos', 'an_sin', 'sa_cos', 'sa_sin', 'RMSresid'}))
            % var
            aus8.(aus8_names{a}).var.(var_names{c}) = ...
                aus8_raw.(aus8_names{a}).var.(var_names{c})(lat_ind, lon_ind, :);
            
        else
            % var
            aus8.(aus8_names{a}).var.(var_names{c}) = ...
                aus8_raw.(aus8_names{a}).var.(var_names{c})(lat_ind, lon_ind, depth_ind);
            
        end
        % att
        aus8.(aus8_names{a}).att.(var_names{c}) = ...
            aus8_raw.(aus8_names{a}).att.(var_names{c});
        
    end
    
    % coors
    aus8.(aus8_names{a}).var.lon = lon;
    aus8.(aus8_names{a}).var.lat = lat;
    aus8.(aus8_names{a}).var.depth = depth;
    aus8.(aus8_names{a}).var.depth_ann = ...
        aus8_raw.(aus8_names{a}).var.depth_ann;
    aus8.(aus8_names{a}).var.depth_san = ...
        aus8_raw.(aus8_names{a}).var.depth_san;
end


%% save aus8_coast
save('cars_out/aus8', 'aus8')
disp('Saved into aus8_coast !')

