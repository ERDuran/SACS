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
aus8_coor.lat = lat;
aus8_coor.lon = lon;
aus8_coor.depth = depth;

aus8_coor.bottom_depth = NaN([length(lat), length(lon)]);

aus8_temp.mean = aus8.temp.mean;
aus8_salt.mean = aus8.salt.mean;

month_names = {...
    'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
for t = 1 : 12
    aus8_temp.(month_names{t}) = aus8_monthly.temp.(month_names{t});
    aus8_salt.(month_names{t}) = aus8_monthly.salt.(month_names{t});
end

for m = 1 : length(lat)
    for n = 1 : length(lon)
        last_finite = ...
            find(isfinite(squeeze(aus8.temp.mean(m,n,:))), ...
            1, 'last');
        % apply ETOPO5 mask
        if ~isempty(last_finite)
            aus8_coor.bottom_depth(m,n) = depth(last_finite);
            if topog(m,n) > aus8_coor.bottom_depth(m,n)
                land_idx = depth < topog(m,n);
                aus8_temp.mean(m,n,land_idx) = NaN;
                aus8_salt.mean(m,n,land_idx) = NaN;
            end
        end
        
        above_first_NaN = ...
            find(isnan(squeeze(aus8.temp.mean(m,n,:))), ...
            1, 'first') -1;
        % interpolate cavities
        if above_first_NaN ~= last_finite
            aus8_temp.mean(m,n,:) = ...
                fixgaps(squeeze(aus8.temp.mean(m,n,:)));
            aus8_salt.mean(m,n,:) = ...
                fixgaps(squeeze(aus8.salt.mean(m,n,:)));
        end
        
        for t = 1 : 12
            if ~isempty(last_finite)
                if topog(m,n) > aus8_coor.bottom_depth(m,n)
                    aus8_temp.(month_names{t})(m,n,land_idx) = NaN;
                    aus8_salt.(month_names{t})(m,n,land_idx) = NaN;
                end
            end
            if above_first_NaN ~= last_finite
                aus8_temp.(month_names{t})(m,n,:) = ...
                    fixgaps(squeeze(...
                    aus8_monthly.temp.(month_names{t})(m,n,:)));
                aus8_salt.(month_names{t})(m,n,:) = ...
                    fixgaps(squeeze(...
                    aus8_monthly.salt.(month_names{t})(m,n,:)));
            end
        end
    end
    fprintf('lat = %5.3f\n', lat(m))
end


%% save updated a8c
save([data_path 'SACS_data/aus8_coor'], 'aus8_coor')
save([data_path 'SACS_data/aus8_temp'], 'aus8_temp')
save([data_path 'SACS_data/aus8_salt'], 'aus8_salt')
disp('aus8_coor aus8_temp aus8_salt DONE')

