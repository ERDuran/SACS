%% Calculate conservative temperature Theta, absolute salinity asal,
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_z_hr_temp'])
load([data_path 'SACS_data/aus8_z_hr_salt'])
load([data_path 'SACS_data/ETOPO5'])

lat = aus8_z_hr_temp.lat;
lon = aus8_z_hr_temp.lon;
depth = aus8_z_hr_temp.depth;
pres = gsw_p_from_z(depth,-40);

aus8_z_hr_thet.lat = lat;
aus8_z_hr_thet.lon = lon;
aus8_z_hr_thet.depth = depth;


%% Calculate conservative temperature, absolute salinity
% make templates
[aus8_z_hr_asal.mean, aus8_z_hr_thet.mean, aus8_z_hr_rho.mean] = ...
    deal(NaN([length(lat) length(lon) length(depth)]));

month_names = {...
    'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
for t = 1 : 12
    aus8_z_hr_asal.(month_names{t}) = ...
        NaN([length(lat), length(lon), length(depth)]);
    aus8_z_hr_thet.(month_names{t}) = ...
        NaN([length(lat), length(lon), length(depth)]);
    aus8_z_hr_rho.(month_names{t}) = ...
        NaN([length(lat), length(lon), length(depth)]);
end

% Calculate using functions from Gibbs seawater toolbox TEOS-10
for p = 1 : length(pres)
    % Mean absolute salinity and conservative temperature from
    % practical salinity and in-situ temperature
    aus8_z_hr_asal.mean(:,:,p) = gsw_SA_from_SP(...
        aus8_z_hr_salt.mean(:,:,p),pres(p),lon,lat);
    aus8_z_hr_thet.mean(:,:,p) = gsw_CT_from_t(...
        aus8_z_hr_asal.mean(:,:,p),...
        aus8_z_hr_temp.mean(:,:,p),pres(p));
    
    for t = 1 : 12
        aus8_z_hr_asal.(month_names{t})(:,:,p) = gsw_SA_from_SP(...
            aus8_z_hr_salt.(month_names{t})(:,:,p),pres(p),lon,lat);
        
        aus8_z_hr_thet.(month_names{t})(:,:,p) = gsw_CT_from_t(...
            aus8_z_hr_asal.(month_names{t})(:,:,p),...
            aus8_z_hr_temp.(month_names{t})(:,:,p),pres(p));
    end
    fprintf('p = %4.0f \n',depth(p))
end


% Calculate density
for p = 1 : length(pres)
    aus8_z_hr_rho.mean(:,:,p) = gsw_rho(...
        aus8_z_hr_asal.mean(:,:,p),...
        aus8_z_hr_thet.mean(:,:,p),pres(p));
    
    for t = 1 : 12
        aus8_z_hr_rho.(month_names{t})(:,:,p) = gsw_rho(...
            aus8_z_hr_asal.(month_names{t})(:,:,p),...
            aus8_z_hr_thet.(month_names{t})(:,:,p),pres(p));
    end
end


%% save the structure
save([data_path 'SACS_data/aus8_z_hr_asal'], 'aus8_z_hr_asal')
save([data_path 'SACS_data/aus8_z_hr_thet'], 'aus8_z_hr_thet')
save([data_path 'SACS_data/aus8_z_hr_rho'], 'aus8_z_hr_rho')

disp('aus8_z_hr_asal aus8_z_hr_thet aus8_z_hr_rho DONE')

