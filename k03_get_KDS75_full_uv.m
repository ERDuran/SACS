%%
clearvars('-except', '*_path')

load([data_path 'SACS_data/aus8_coor'])
load([data_path 'SACS_data/SmSan02'])


MTH = {'JFM', 'AMJ', 'JAS', 'OND'};
aus8_coor.MTH = MTH;
save([data_path 'SACS_data/aus8_coor'], 'aus8_coor')
lat = aus8_coor.lat;
lon = aus8_coor.lon;
lat_u = aus8_coor.lat_u;
lon_u = aus8_coor.lon_u;
lat_v = aus8_coor.lat_v;
lon_v = aus8_coor.lon_v;
depth_mid = aus8_coor.depth_mid;
depth = aus8_coor.depth;


%%
for t = 1 : 4
    [data.(MTH{t}), ~, kds_att] = ...
        nc2mat([data_path ...
        'KDS75/KDS75_ncra_y103to109_' MTH{t} '.nc'], ...
        'ALL');
        
    disp([MTH{t} ' OK!'])
end
% fix axis
yt_ocean = flipud(data.JFM.yt_ocean);
xt_ocean = (data.JFM.xt_ocean + 360)';
st_ocean = -data.JFM.st_ocean;
yu_ocean = flipud(data.JFM.yu_ocean);
xu_ocean = (data.JFM.xu_ocean + 360)';
sw_ocean = -data.JFM.sw_ocean;

% u, st_ocean, yu_ocean, xu_ocean
% v, st_ocean, yu_ocean, xu_ocean
% wt, sw_ocean, yt_ocean, xt_ocean
% tx_trans, st_ocean, yt_ocean, xu_ocean
% ty_trans, st_ocean, yu_ocean, xt_ocean



%%
KDS75_fulu.mean = 0;
KDS75_fulv.mean = 0;
KDS75_fulw.mean = 0;

for t = 1 : 4
    KDS75_fulu.(MTH{t}) = ...
        flipud(permute(data.(MTH{t}).u, [2, 1, 3]));
    KDS75_fulv.(MTH{t}) = ...
        flipud(permute(data.(MTH{t}).v, [2, 1, 3]));
    KDS75_fulw.(MTH{t}) = ...
        flipud(permute(data.(MTH{t}).wt, [2, 1, 3]));
    
    KDS75_fulu.(MTH{t})(isnan(KDS75_fulu.(MTH{t}))) = 0;
    KDS75_fulv.(MTH{t})(isnan(KDS75_fulv.(MTH{t}))) = 0;
    KDS75_fulw.(MTH{t})(isnan(KDS75_fulw.(MTH{t}))) = 0;
    
end
mean_u = KDS75_fulu.JFM;
mean_v = KDS75_fulv.JFM;
mean_w = KDS75_fulw.JFM;
for t = 2 : 4
    mean_u = mean_u + KDS75_fulu.(MTH{t});
    mean_v = mean_v + KDS75_fulv.(MTH{t});
    mean_w = mean_w + KDS75_fulw.(MTH{t});
end
KDS75_fulu.mean = mean_u/4;
KDS75_fulv.mean = mean_v/4;
KDS75_fulw.mean = mean_w/4;


%% interp2: same lat/lon first
KDau_fulu_dum.mean = ...
    NaN([length(lat_u) length(lon_u) length(st_ocean)]);
KDau_fulv_dum.mean = ...
    NaN([length(lat_v) length(lon_v) length(st_ocean)]);
KDau_fulw_dum.mean = ...
    deal(NaN([length(lat_u) length(lon_v) length(sw_ocean)]));

for t = 1 : 4
    KDau_fulu_dum.(MTH{t}) =  ...
        NaN([length(lat_u), length(lon_u), length(st_ocean)]);
    KDau_fulv_dum.(MTH{t}) =  ...
        NaN([length(lat_v), length(lon_v), length(st_ocean)]);
    KDau_fulw_dum.(MTH{t}) = ...
        NaN([length(lat_u), length(lon_v), length(sw_ocean)]);
end

for p = 1 : length(st_ocean)
    KDau_fulu_dum.mean(:,:,p) = interp2(xu_ocean, yu_ocean, ...
        KDS75_fulu.mean(:,:,p), lon_u, lat_u);
    KDau_fulv_dum.mean(:,:,p) = interp2(xu_ocean, yu_ocean, ...
        KDS75_fulv.mean(:,:,p), lon_v, lat_v);
    KDau_fulw_dum.mean(:,:,p) = interp2(xt_ocean, yt_ocean, ...
        KDS75_fulw.mean(:,:,p), lon_v, lat_u);
    
    for t = 1 : 4
        KDau_fulu_dum.(MTH{t})(:,:,p) = interp2(xu_ocean, yu_ocean, ...
            KDS75_fulu.(MTH{t})(:,:,p), lon_u, lat_u);
        KDau_fulv_dum.(MTH{t})(:,:,p) = interp2(xu_ocean, yu_ocean, ...
            KDS75_fulv.(MTH{t})(:,:,p), lon_v, lat_v);
        KDau_fulw_dum.(MTH{t})(:,:,p) = interp2(xt_ocean, yt_ocean, ...
            KDS75_fulw.(MTH{t})(:,:,p), lon_v, lat_u);
    end
    fprintf('depth = %7.3f \n', st_ocean(p))
end


%% interp1: same depth second (this takes a long time to run!)
% make templates
KDau_fulu.mean = ...
    NaN([length(lat_u) length(lon_u) length(depth_mid)]);
KDau_fulv.mean = ...
    NaN([length(lat_v) length(lon_v) length(depth_mid)]);
KDau_fulw.mean = NaN([length(lat_u) length(lon_v) length(depth)]);

for t = 1 : 4
    KDau_fulu.(MTH{t}) = ...
        NaN([length(lat_u), length(lon_u), length(depth_mid)]);
    KDau_fulv.(MTH{t}) = ...
        NaN([length(lat_v), length(lon_v), length(depth_mid)]);
    KDau_fulw.(MTH{t}) = NaN([length(lat_u) length(lon_v) length(depth)]);
end

for m = 1 : length(lat_u)
    for n = 1 : length(lon_u)
        KDau_fulu.mean(m,n,:) = interp1(st_ocean, ...
            squeeze(KDau_fulu_dum.mean(m,n,:)), depth_mid);
        
        for t = 1 : 4
            KDau_fulu.(MTH{t})(m,n,:) = interp1(st_ocean, ...
                squeeze(KDau_fulu_dum.(MTH{t})(m,n,:)), depth_mid);
        end
    end
    fprintf('lat = %7.3f \n',lat(m))
end


for m = 1 : length(lat_v)
    for n = 1 : length(lon_v)
        KDau_fulv.mean(m,n,:) = interp1(st_ocean, ...
            squeeze(KDau_fulv_dum.mean(m,n,:)), depth_mid);
        
        for t = 1 : 4
            KDau_fulv.(MTH{t})(m,n,:) = interp1(st_ocean, ...
                squeeze(KDau_fulv_dum.(MTH{t})(m,n,:)), depth_mid);
        end
    end
    fprintf('lat = %7.3f \n',lat(m))
end


for m = 1 : length(lat_u)
    for n = 1 : length(lon_v)
        KDau_fulw.mean(m,n,:) = interp1(sw_ocean, ...
            squeeze(KDau_fulw_dum.mean(m,n,:)), depth);
        
        for t = 1 : 4
            KDau_fulw.(MTH{t})(m,n,:) = interp1(sw_ocean, ...
                squeeze(KDau_fulw_dum.(MTH{t})(m,n,:)), depth);
        end
    end
    fprintf('lat = %7.3f \n',lat(m))
end


%% apply CARS u,v mask.
MTH{5} = 'mean';
for t = 1 : 5
    KDau_fulu.(MTH{t})(KDau_fulu.(MTH{t}) == 0) = NaN;
    KDau_fulv.(MTH{t})(KDau_fulv.(MTH{t}) == 0) = NaN;
    KDau_fulw.(MTH{t})(KDau_fulw.(MTH{t}) == 0) = NaN;
    
    for m = 1 : length(lat_u)
        for n = 1 : length(lon_u)
            last_finite = ...
                find(isfinite(squeeze(KDau_fulu.(MTH{t})(m,n,:))), ...
                1, 'last');
            if ~isempty(last_finite)
                if depth_mid(last_finite) < aus8_coor.u_bottom(m,n)
                    land_idx = depth_mid < aus8_coor.u_bottom(m,n);
                    KDau_fulu.(MTH{t})(m,n,land_idx) = NaN;
                end
            end
        end
    end
    for m = 1 : length(lat_v)
        for n = 1 : length(lon_v)
            last_finite = ...
                find(isfinite(squeeze(KDau_fulv.(MTH{t})(m,n,:))), ...
                1, 'last');
            if ~isempty(last_finite)
                if depth_mid(last_finite) < aus8_coor.v_bottom(m,n)
                    land_idx = depth_mid < aus8_coor.v_bottom(m,n);
                    KDau_fulv.(MTH{t})(m,n,land_idx) = NaN;
                end
            end
        end
    end
    
    
    
    for m = 1 : length(lat_u)
        for n = 1 : length(lon_v)
            last_finite = ...
                find(isfinite(squeeze(KDau_fulw.(MTH{t})(m,n,:))), ...
                1, 'last');
            if ~isempty(last_finite)
                current_uv = [aus8_coor.u_bottom(m,n:n+1), ...
                    aus8_coor.v_bottom(m:m+1,n)'];
                deepest_uv = min(current_uv);
                
                if deepest_uv == depth(end)
                    continue
                    
                elseif deepest_uv >= depth(last_finite)  
                    land_idx = depth <= deepest_uv;
                    KDau_fulw.(MTH{t})(m,n,land_idx) = NaN;
                    
                end
            end
        end
    end
    disp([MTH{t} ' OK!'])
end


%%
save([data_path 'SACS_data/KDau_fulu'], 'KDau_fulu')
save([data_path 'SACS_data/KDau_fulv'], 'KDau_fulv')
save([data_path 'SACS_data/KDau_fulw'], 'KDau_fulw')
disp('KDau_fulu KDau_fulv KDau_fulw DONE')

