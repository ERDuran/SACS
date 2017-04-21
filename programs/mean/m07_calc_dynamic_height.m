%% Calculate dynamic height and geostrophic velocities in each AOI
clearvars('-except', '*_path')

load aus8_TEOS10
lat = aus8_TEOS10.lat;
lon = aus8_TEOS10.lon;
pres = aus8_TEOS10.pres;
asal = aus8_TEOS10.asal.mean;
Theta = aus8_TEOS10.Theta.mean;


%% Calculate dynamic height
disp('Calculating dynamic heights at the surface...')

% pressure reference: the depth of no motion
p_ref = 0;

%
aus8_streamfunction.lat = lat;
aus8_streamfunction.lon = lon;
aus8_streamfunction.pres = pres;
aus8_streamfunction.p_ref_list = p_ref;
aus8_streamfunction.dynh.name = 'Dynamic Height';
aus8_streamfunction.dynh.symbol = 'psi';
aus8_streamfunction.dynh.units = 'm^2/s^2';

for pp = 1 : length(p_ref)
    % template
    dynh = NaN(size(Theta));

    for ii = 1 : length(lat)
        for jj = 1 : length(lon)
            
            % find first NaN point in the current water column
            first_NaN_ind = find(isnan(Theta(ii,jj,:)), 1, 'first');
            
            % if the land is shallower or at the pressure reference
            if pres(first_NaN_ind) <= p_ref(pp)
                % then we can't calculate the streamfunction there
                continue
                
                % otherwise, we're good
            else
                % calculate the streamfunction at the surface
                dynh(ii,jj,:) = gsw_geo_strf_dyn_height(...
                    squeeze(asal(ii,jj,:)), ...
                    squeeze(Theta(ii,jj,:)), ...
                    pres, p_ref(pp));
            end
        end
        
        fprintf('dynh with p_ref = %3.0f at lat = %7.3f OK ! \n', ...
            p_ref(pp), lat(ii))
    end
    
    %
    aus8_streamfunction.dynh.(['p_ref_' num2str(p_ref(pp))]) = dynh;
    fprintf('Dynamic Height Calculated OK \n')
end


%%
save([cars_out_path 'aus8_streamfunction'], 'aus8_streamfunction')
disp(['aus8_streamfunction saved in ' ...
    cars_out_path 'aus8_streamfunction'])

