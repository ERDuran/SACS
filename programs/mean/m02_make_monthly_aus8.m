%% make monthly variations
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8'])


%% 
time_step = [...
    15, 46, 74, 105, 135, 166, ...
    196, 227, 258, 288, 319, 349];

month_names = {...
    'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

% time cycle setup
time_serie = 2 * pi / 365 * time_step;

for ll = 1 : 12
    temp_monthly = NaN(size(aus8.temp.mean));
    salt_monthly = NaN(size(aus8.salt.mean));
    
    for kk = 1 : length(aus8.depth)
        if aus8.depth(kk) >= aus8.depth_semiann(end)
            temp_monthly(:,:,kk) = ...
                ... % mean component
                aus8.temp.mean(:,:,kk) + ...
                ... % annual component
                aus8.temp.an_cos(:,:,kk)*cos(time_serie(ll)) + ...
                aus8.temp.an_sin(:,:,kk)*sin(time_serie(ll)) + ...
                ... % semi-annual component
                aus8.temp.sa_cos(:,:,kk)*cos(2*time_serie(ll)) + ...
                aus8.temp.sa_sin(:,:,kk)*sin(2*time_serie(ll));
            
            salt_monthly(:,:,kk) = ...
                aus8.salt.mean(:,:,kk) + ...
                aus8.salt.an_cos(:,:,kk)*cos(time_serie(ll)) + ...
                aus8.salt.an_sin(:,:,kk)*sin(time_serie(ll)) + ...
                aus8.salt.sa_cos(:,:,kk)*cos(2*time_serie(ll)) + ...
                aus8.salt.sa_sin(:,:,kk)*sin(2*time_serie(ll));
            
            % depths above the annual component max depth
        elseif aus8.depth(kk) >= aus8.depth_ann(end)
            temp_monthly(:,:,kk) = ...
                ... % mean component
                aus8.temp.mean(:,:,kk) + ...
                ... % annual component
                aus8.temp.an_cos(:,:,kk)*cos(time_serie(ll)) + ...
                aus8.temp.an_sin(:,:,kk)*sin(time_serie(ll));
            
            salt_monthly(:,:,kk) = ...
                aus8.salt.mean(:,:,kk) + ...
                aus8.salt.an_cos(:,:,kk)*cos(time_serie(ll)) + ...
                aus8.salt.an_sin(:,:,kk)*sin(time_serie(ll));
            
            % anything below the max depth of the annual component
            % is the simply the mean
        else
            temp_monthly(:,:,kk) = ...
                ... % mean component
                aus8.temp.mean(:,:,kk);
            
            salt_monthly(:,:,kk) = ...
                aus8.salt.mean(:,:,kk);
        end
    end
    
    aus8_monthly.temp.(month_names{ll}) = temp_monthly;
    aus8_monthly.salt.(month_names{ll}) = salt_monthly;
    disp([month_names{ll} ' OK!'])
end


%% save aus8_coast
save([data_path 'SACS_data/aus8_monthly'], 'aus8_monthly')
disp('aus8_monthly DONE')

