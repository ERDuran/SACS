%% calc t s and rho
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8'])
load([data_path 'SACS_data/aus8_monthly'])
load([data_path 'SACS_data/ETOPO5'])

lat = aus8.lat;
lon = aus8.lon;
depth = aus8.depth;
topog = ETOPO5.topo_sm_interp;


%% increase vertical resolution
depth_hr = 0 : -10 : aus8.depth(end);
aus8_z_hr_temp.lat = lat;
aus8_z_hr_temp.lon = lon;
aus8_z_hr_temp.depth = depth_hr;

aus8_z_hr_temp.bottom_depth = NaN([length(lat), length(lon)]);

aus8_z_hr_temp.mean = ...
    NaN([length(lat), length(lon), length(depth_hr)]);
aus8_z_hr_salt.mean = ...
    NaN([length(lat), length(lon), length(depth_hr)]);

month_names = {...
    'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
for t = 1 : 12
    aus8_z_hr_temp.(month_names{t}) = ...
        NaN([length(lat), length(lon), length(depth_hr)]);
    aus8_z_hr_salt.(month_names{t}) = ...
        NaN([length(lat), length(lon), length(depth_hr)]);
end

for m = 1 : length(lat)
    for n = 1 : length(lon)
        aus8_z_hr_temp.mean(m,n,:) = ...
            interp1(depth,squeeze(aus8.temp.mean(m,n,:)),depth_hr);
        aus8_z_hr_salt.mean(m,n,:) = ...
            interp1(depth,squeeze(aus8.temp.mean(m,n,:)),depth_hr);
        
        % interpolate cavities
        above_first_NaN = ...
            find(isnan(squeeze(aus8_z_hr_temp.mean(m,n,:))), ...
            1, 'first') -1;
        last_finite = ...
            find(isfinite(squeeze(aus8_z_hr_temp.mean(m,n,:))), ...
            1, 'last');
        if above_first_NaN ~= last_finite
            aus8_z_hr_temp.mean(m,n,:) = ...
                fixgaps(squeeze(aus8_z_hr_temp.mean(m,n,:)));
            aus8_z_hr_salt.mean(m,n,:) = ...
                fixgaps(squeeze(aus8_z_hr_salt.mean(m,n,:)));
        end
        
        % apply ETOPO5 mask
        if ~isempty(last_finite)
            aus8_z_hr_temp.bottom_depth(m,n) = depth_hr(last_finite);
            if topog(m,n) > aus8_z_hr_temp.bottom_depth(m,n)
                land_idx = depth_hr < topog(m,n);
                aus8_z_hr_temp.mean(m,n,land_idx) = NaN;
                aus8_z_hr_salt.mean(m,n,land_idx) = NaN;
            end
        end
        
        for t = 1 : 12
            aus8_z_hr_temp.(month_names{t})(m,n,:) = ...
                interp1(depth,squeeze(aus8.temp.mean(m,n,:)),depth_hr);
            aus8_z_hr_salt.(month_names{t})(m,n,:) = ...
                interp1(depth,squeeze(aus8.temp.mean(m,n,:)),depth_hr);
            
            if above_first_NaN ~= last_finite
                aus8_z_hr_temp.(month_names{t})(m,n,:) = ...
                    fixgaps(squeeze(...
                    aus8_z_hr_temp.(month_names{t})(m,n,:)));
                aus8_z_hr_salt.(month_names{t})(m,n,:) = ...
                    fixgaps(squeeze(...
                    aus8_z_hr_salt.(month_names{t})(m,n,:)));
            end
            
            if ~isempty(last_finite)
                if topog(m,n) > aus8_z_hr_temp.bottom_depth(m,n)
                    aus8_z_hr_temp.(month_names{t})(m,n,land_idx) = NaN;
                    aus8_z_hr_salt.(month_names{t})(m,n,land_idx) = NaN;
                end
            end
        end
    end
    fprintf('lat = %5.3f\n', lat(m))
end


%% save updated a8c
save([data_path 'SACS_data/aus8_z_hr_temp'], 'aus8_z_hr_temp')
save([data_path 'SACS_data/aus8_z_hr_salt'], 'aus8_z_hr_salt')
disp('aus8_z_hr_temp aus8_z_hr_salt DONE')

