%% Calculate conservative temperature Theta, absolute salinity asal,
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_z_hr_temp'])
load([data_path 'SACS_data/aus8_z_hr_salt'])
load([data_path 'SACS_data/ETOPO5'])

lat = aus8_z_hr_temp.lat;
lon = aus8_z_hr_temp.lon;
depth = aus8_z_hr_temp.depth;
pres = gsw_p_from_z(depth,-40);

aus8_z_hr_Theta.lat = lat;
aus8_z_hr_Theta.lon = lon;
aus8_z_hr_Theta.depth = depth;


%% Calculate conservative temperature, absolute salinity
disp('Calculating Theta and asal ...')

% make templates
[asal, Theta] = deal(NaN([length(lat) length(lon) length(pres)]));

% Calculate using functions from Gibbs seawater toolbox TEOS-10
for kk = 1 : length(pres)
    % Mean absolute salinity and conservative temperature from
    % practical salinity and in-situ temperature
    asal(:,:,kk) = gsw_SA_from_SP(sal(:,:,kk),pres(kk),lon,lat);
    Theta(:,:,kk) = gsw_CT_from_t(asal(:,:,kk),temp(:,:,kk),pres(kk));
    
    fprintf('p = %4.0f \n',pres(kk))
end
clear temp sal page


% interpolate lateral bumps above "cavities: Theta and asal
% ie. "cavities" T and S: ocean points that are underneath land points
Theta_fixed = Theta;
asal_fixed = asal;

cav = NaN(length(lat), length(lon));

for ii = 1 : length(lat)
    for jj = 1 : length(lon)
        above_first_NaN = ...
            find(isnan(squeeze(Theta(ii,jj,:))), 1, 'first') -1;
        last_finite = ...
            find(isfinite(squeeze(Theta(ii,jj,:))), 1, 'last');
        
%         if isempty(above_first_NaN) || isempty(last_finite)
%             continue
            
        if above_first_NaN ~= last_finite
            cav(ii,jj) = 1;
            Theta_fixed(ii,jj,:) = fixgaps(squeeze(Theta(ii,jj,:)));
            asal_fixed(ii,jj,:) = fixgaps(squeeze(asal(ii,jj,:)));
        end
        
    end
end


% mask data with ETOPO-5
topog = topog_mask.aus8.topog;

asal_masked = asal_fixed;
Theta_masked = Theta_fixed;

% apply mask on temperature and salinity
for ii = 1 : length(lat)
    for jj = 1 : length(lon)
        last_finite = find(isfinite(Theta_masked(ii,jj,:)),1,'last');
        
        if ~isempty(last_finite)
            topog_now = topog(ii,jj);
            deepest_now = -pres(last_finite);
            
            if topog_now > deepest_now
                land_idx = -pres < topog_now;
                asal_masked(ii,jj,land_idx) = NaN;
                Theta_masked(ii,jj,land_idx) = NaN;
                
            end
        end
    end
end


aus8_TEOS10.asal.name = 'Absolute Salinity';
aus8_TEOS10.asal.symbol = 'S_A';
aus8_TEOS10.asal.units = 'g/kg';
aus8_TEOS10.asal.mean = asal_masked;
aus8_TEOS10.Theta.name = 'Conservative Temperature';
aus8_TEOS10.Theta.symbol = '\Theta';
aus8_TEOS10.Theta.units = '\circ C';
aus8_TEOS10.Theta.mean = Theta_masked;
disp('Theta and asal calculated and saved into cars_out')


%% Calculate potential density anomaly
disp('Calculating rho ...')

rho = NaN([length(lat) length(lon) length(pres)]);

% in situ density
for kk = 1 : length(pres)
    rho(:,:,kk) = gsw_rho(asal_fixed(:,:,kk),Theta_fixed(:,:,kk),pres(kk));
end

aus8_TEOS10.rho.name = 'Density';
aus8_TEOS10.rho.symbol = '\rho';
aus8_TEOS10.rho.units = 'kg/m^3';
aus8_TEOS10.rho.mean = rho;

disp('rho calculated and saved into cars_out')


%% save the structure
save('cars_out/aus8_TEOS10', 'aus8_TEOS10')
disp('Saved in aus8_TEOS10.')

