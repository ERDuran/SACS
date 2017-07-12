%%
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/KDau_thet'])
load([data_path 'SACS_data/KDau_salt'])

lat = aus8_coor.lat;
lon = aus8_coor.lon;
depth = aus8_coor.depth;
Months = aus8_coor.Months;
pres = aus8_coor.pres;


%% Calculate conservative temperature, absolute salinity
% make templates
[KDau_asal.mean, KDau_rho.mean] = ...
    deal(NaN([length(lat) length(lon) length(depth)]));
for t = 1 : 12
    KDau_asal.(Months{t}) = ...
        NaN([length(lat), length(lon), length(depth)]);
    KDau_rho.(Months{t}) = ...
        NaN([length(lat), length(lon), length(depth)]);
end

% Calculate using functions from Gibbs seawater toolbox TEOS-10
for p = 1 : length(pres)
    % Mean absolute salinity and conservative temperature from
    % practical salinity and in-situ temperature
    KDau_asal.mean(:,:,p) = gsw_SA_from_SP(...
        KDau_salt.mean(:,:,p),pres(p),lon,lat);
    
    for t = 1 : 12
        KDau_asal.(Months{t})(:,:,p) = gsw_SA_from_SP(...
            KDau_salt.(Months{t})(:,:,p),pres(p),lon,lat);
    end
    fprintf('z = %4.0f \n',depth(p))
end


% Calculate density
for p = 1 : length(pres)
    KDau_rho.mean(:,:,p) = gsw_rho(...
        KDau_asal.mean(:,:,p),...
        KDau_thet.mean(:,:,p),pres(p));
    
    for t = 1 : 12
        KDau_rho.(Months{t})(:,:,p) = gsw_rho(...
            KDau_asal.(Months{t})(:,:,p),...
            KDau_thet.(Months{t})(:,:,p),pres(p));
    end
end


%%
save([data_path 'SACS_data/KDau_asal'], 'KDau_asal')
save([data_path 'SACS_data/KDau_rho'], 'KDau_rho')
disp('KDau_asal KDau_rho DONE')

