%% Calculate dynamic height and geostrophic velocities in each AOI
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/aus8_thet'])
load([data_path 'SACS_data/aus8_asal'])

lat = aus8_coor.lat;
lon = aus8_coor.lon;
depth = aus8_coor.depth;
pres = aus8_coor.pres;


%% Calculate dynamic height
aus8_dynh0.mean = NaN(size(aus8_thet.mean));
month_names = {...
    'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
for t = 1 : 12
    aus8_dynh0.(month_names{t}) = ...
        NaN(size(aus8_thet.mean));
end

for m = 1 : length(lat)
    for n = 1 : length(lon)
        
        % find first NaN point in the current water column
        first_NaN_ind = ...
            find(isnan(aus8_thet.mean(m,n,:)), 1, 'first');
        
        % if the land is shallower or at the pressure reference
        if depth(first_NaN_ind) >= 0
            % then we can't calculate the streamfunction there
            continue
            
            % otherwise, we're good
        else
            % calculate the streamfunction at the surface
            aus8_dynh0.mean(m,n,:) = gsw_geo_strf_dyn_height(...
                squeeze(aus8_asal.mean(m,n,:)), ...
                squeeze(aus8_thet.mean(m,n,:)), ...
                pres, 0);
            
            for t = 1 : 12
                aus8_dynh0.(month_names{t})(m,n,:) = ...
                    gsw_geo_strf_dyn_height(...
                    squeeze(aus8_asal.(month_names{t})(m,n,:)), ...
                    squeeze(aus8_thet.(month_names{t})(m,n,:)), ...
                    pres, 0);
            end
        end
    end
    
    fprintf('dynh at lat = %7.3f OK ! \n', ...
        lat(m))
end


%% smooth out dynamic height
for m = 2 : length(lat)-1
    for n = 2 : length(lon)-1
        finite_dynh = find(isfinite(aus8_dynh0.mean(m,n,:)));
        
        if ~isempty(finite_dynh)
            %
            dynh_raw_N = squeeze(aus8_dynh0.mean(m-1,n,finite_dynh));
            dynh_raw_S = squeeze(aus8_dynh0.mean(m+1,n,finite_dynh));
            dynh_raw_W = squeeze(aus8_dynh0.mean(m,n-1,finite_dynh));
            dynh_raw_E = squeeze(aus8_dynh0.mean(m,n+1,finite_dynh));
            dynh_raw_C = squeeze(aus8_dynh0.mean(m,n,finite_dynh));
            
            %
            % %%% 1 way: if any of the four points around are nan then
            % the centre point is the original value ie. only averaging
            % over the 5 points
            % PRETTY GOOD WAY
            %
            % mean_dynh = mean([...
            %     dynh_raw_N, ...
            %     dynh_raw_S, ...
            %     dynh_raw_W, ...
            %     dynh_raw_E, ...
            %     dynh_raw_C],2);
            %
            % mean_dynh_NaN = isnan(mean_dynh);
            % mean_dynh(mean_dynh_NaN) = dynh_raw_C(mean_dynh_NaN);
            % dynh(ii,jj,finite_dynh) = mean_dynh;
            % %%%
            %
            %
            % %%% 2 way: average along the whole finite water column of
            % % the centre ie. allowed to average over less than 5
            % % points BAD WAY
            % mean_dynh = nanmean([...
            %     dynh_raw_N, ...
            %     dynh_raw_S, ...
            %     dynh_raw_W, ...
            %     dynh_raw_E, ...
            %     dynh_raw_C],2);
            %
            % %
            % dynh(ii,jj,finite_dynh) = mean_dynh;
            % %%%
            
            
            %%% 3 way: if we are above the shelf (ie. depth shallower
            % than 2000 db), then do not average at all and take
            % original centre cell. the average offshore is always over
            % 5 points
            % SLIGHTLY BETTER THAN 1
            if finite_dynh(end) == length(pres)
                %
                mean_dynh = mean([...
                    dynh_raw_N, ...
                    dynh_raw_S, ...
                    dynh_raw_W, ...
                    dynh_raw_E, ...
                    dynh_raw_C],2);
                %
                mean_dynh_NaN = isnan(mean_dynh);
                mean_dynh(mean_dynh_NaN) = dynh_raw_C(mean_dynh_NaN);
                aus8_dynh0.mean(m,n,finite_dynh) = mean_dynh;
            else
                %
                aus8_dynh0.mean(m,n,finite_dynh) = dynh_raw_C;
            end
            %%%
            
            %
            % %%% 4 way: if we are above the shelf (ie. depth shallower
            % % than 2000 db), then do not average at all and take
            % % original centre cell. the average offshore can be over
            % % less than 5 points
            % % PRETTY BAD WAY
            % if finite_dynh(end) == length(pres)
            %     %
            %     mean_dynh = nanmean([...
            %         dynh_raw_N, ...
            %         dynh_raw_S, ...
            %         dynh_raw_W, ...
            %         dynh_raw_E, ...
            %         dynh_raw_C],2);
            %
            %     %
            %     dynh(ii,jj,finite_dynh) = mean_dynh;
            %
            % else
            %     %
            %     dynh(ii,jj,finite_dynh) = dynh_raw_C;
            %
            % end
            % %%%
            
            
            % monthly data
            for t = 1 : 12
                dynh_raw_N = squeeze(...
                    aus8_dynh0.(month_names{t})(m-1,n,finite_dynh));
                dynh_raw_S = squeeze(...
                    aus8_dynh0.(month_names{t})(m+1,n,finite_dynh));
                dynh_raw_W = squeeze(...
                    aus8_dynh0.(month_names{t})(m,n-1,finite_dynh));
                dynh_raw_E = squeeze(...
                    aus8_dynh0.(month_names{t})(m,n+1,finite_dynh));
                dynh_raw_C = squeeze(...
                    aus8_dynh0.(month_names{t})(m,n,finite_dynh));
                
                if finite_dynh(end) == length(pres)
                    mean_dynh = mean([...
                        dynh_raw_N, ...
                        dynh_raw_S, ...
                        dynh_raw_W, ...
                        dynh_raw_E, ...
                        dynh_raw_C],2);
                    mean_dynh_NaN = isnan(mean_dynh);
                    mean_dynh(mean_dynh_NaN) = dynh_raw_C(mean_dynh_NaN);
                    aus8_dynh0.(month_names{t})(m,n,finite_dynh) = ...
                        mean_dynh;
                else
                    aus8_dynh0.(month_names{t})(m,n,finite_dynh) = ...
                        dynh_raw_C;
                end
            end
        end
    end
end


%%
save([data_path 'SACS_data/aus8_dynh0'], 'aus8_dynh0')

disp('aus8_dynh0 DONE')

