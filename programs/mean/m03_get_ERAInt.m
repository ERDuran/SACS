%% get wind data
clearvars('-except', '*_path')

% wind
[wind_data, wind_gen_att, wind_att] = ...
    nc2mat([data_path ...
    'ERA_INTERIM/mth_av_turb_surf_stress_20032012/' ...
    '_grib2netcdf-atls12-a562cefde8a29a7288fa0b8b7f9413f7-qJ0lae.nc'], ...
    'ALL');


%%
ERAInt.lat = double(wind_data.latitude);
ERAInt.lon = double(wind_data.longitude)';
ERAInt.time_num = double(wind_data.time);
ERAInt.time_vec = datestr(datenum(1900,1,1) + ERAInt.time_num/24);

%
inss = permute(wind_data.inss, [2,1,3]);
ERAInt.Tau_y.mean = mean(inss,3);

%
iews = permute(wind_data.iews, [2,1,3]);
ERAInt.Tau_x.mean = mean(iews,3);

%
month_names = {...
    'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
month_vector = 1 : 12 : length(ERAInt.time_num);

for ll = 1 : 12
    ERAInt.Tau_y.(month_names{ll}) = ...
        mean(inss(:,:,month_vector+ll-1),3);
    ERAInt.Tau_x.(month_names{ll}) = ...
        mean(iews(:,:,month_vector+ll-1),3);
end


%% save
save([data_path 'SACS_data/ERAInt'], 'ERAInt')
disp('ERAInt DONE')

