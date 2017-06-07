%% Calculate conservative temperature Theta, absolute salinity asal,
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_temp'])
load([data_path 'SACS_data/aus8_salt'])
load([data_path 'SACS_data/ETOPO5'])

lat = aus8_coor.lat;
lon = aus8_coor.lon;
depth = aus8_coor.depth;
pres = gsw_p_from_z(depth,-40);
aus8_coor.pres = pres;
save([data_path 'SACS_data/aus8_coor'], 'aus8_coor')


%% Calculate conservative temperature, absolute salinity
% make templates
[aus8_asal.mean, aus8_thet.mean, aus8_rho.mean] = ...
    deal(NaN([length(lat) length(lon) length(depth)]));

month_names = {...
    'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
for t = 1 : 12
    aus8_asal.(month_names{t}) = ...
        NaN([length(lat), length(lon), length(depth)]);
    aus8_thet.(month_names{t}) = ...
        NaN([length(lat), length(lon), length(depth)]);
    aus8_rho.(month_names{t}) = ...
        NaN([length(lat), length(lon), length(depth)]);
end

% Calculate using functions from Gibbs seawater toolbox TEOS-10
for p = 1 : length(pres)
    % Mean absolute salinity and conservative temperature from
    % practical salinity and in-situ temperature
    aus8_asal.mean(:,:,p) = gsw_SA_from_SP(...
        aus8_salt.mean(:,:,p),pres(p),lon,lat);
    aus8_thet.mean(:,:,p) = gsw_CT_from_t(...
        aus8_asal.mean(:,:,p),...
        aus8_temp.mean(:,:,p),pres(p));
    
    for t = 1 : 12
        aus8_asal.(month_names{t})(:,:,p) = gsw_SA_from_SP(...
            aus8_salt.(month_names{t})(:,:,p),pres(p),lon,lat);
        
        aus8_thet.(month_names{t})(:,:,p) = gsw_CT_from_t(...
            aus8_asal.(month_names{t})(:,:,p),...
            aus8_temp.(month_names{t})(:,:,p),pres(p));
    end
    fprintf('p = %4.0f \n',depth(p))
end


% Calculate density
for p = 1 : length(pres)
    aus8_rho.mean(:,:,p) = gsw_rho(...
        aus8_asal.mean(:,:,p),...
        aus8_thet.mean(:,:,p),pres(p));
    
    for t = 1 : 12
        aus8_rho.(month_names{t})(:,:,p) = gsw_rho(...
            aus8_asal.(month_names{t})(:,:,p),...
            aus8_thet.(month_names{t})(:,:,p),pres(p));
    end
end


%% save the structure
save([data_path 'SACS_data/aus8_asal'], 'aus8_asal')
save([data_path 'SACS_data/aus8_thet'], 'aus8_thet')
save([data_path 'SACS_data/aus8_rho'], 'aus8_rho')

disp('aus8_asal aus8_thet aus8_rho DONE')

