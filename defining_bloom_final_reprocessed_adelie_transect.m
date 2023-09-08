%% Determining Centers and Edges of Blooms from ACROBAT data
% copied from 'defining_bloom_final_reprocessed_linebyline.m'
% repurposed to create threshold for patch ID from entire day of ACROBAT
% survey
% 2Nov22 major edits to optomize & clean edge definition
% 10July2023 added mlparticle to patchID struct for continuous comparison
% of ml particle and FTLE (see 'continuous_comparison.m')
% 28Aug2023 added march 6th and march 10th surveys

% plot results of binary patch and edge definition in 'plot_patch_definition_daily.m'

% data struct created by 'Acrobat_create_struct.m'
load '/Users/jveatch/Documents/MATLAB/SWARM/ACROBAT/data_reprocessed/ACRO_reprocessed.mat'
load '/Volumes/T7_Shield/jmv208/SWARM_data/SWARM_CODAR.mat'
addpath '/Users/jveatch/Documents/MATLAB/SWARM/ACROBAT/CODE' % location of dist_lat_lon

survey_days = ['15-Jan-2020'; '18-Jan-2020'; '21-Jan-2020'; '24-Jan-2020'; 
    '28-Jan-2020'; '01-Feb-2020';'05-Feb-2020';'07-Feb-2020';
    '12-Feb-2020'; '14-Feb-2020'; '18-Feb-2020'; '21-Feb-2020';'22-Feb-2020'; 
    '25-Feb-2020'; '28-Feb-2020'; '03-Mar-2020'; '06-Mar-2020'; '10-Mar-2020' ];
% took out March 6th becuase there is no HFR data there

% load 'sml_max_part.mat';
%% Small detour to calculate density difference
for i = 1:length(ACRO.mld)
    if isnan(ACRO.mld(i))
        ACRO.dens_diff(i) = NaN;
    else
        ind_no_nan = ~isnan(ACRO.dens(:,i));
        dens = ACRO.dens(:,i);
        dens_no_nan = dens(ind_no_nan);
        
        if dens_no_nan(end) < 20 % if the profile did not go below 20 meters
            ACRO.dens_diff(i) = NaN;
        end
        
        ACRO.dens_diff(i) = real(dens_no_nan(end) - dens_no_nan(2));
        
    end
end

%% 
[x,y] = meshgrid(CODAR.lon, CODAR.lat);
codar_line_lat = y(:,1);
codar_line_lon = x(:,1);

dist_codar = nan(size(codar_line_lon));

for i = 1:length(codar_line_lat)-1
    [dist_codar(i), AF, AR] = dist_lat_lon([codar_line_lat(i); codar_line_lat(i+1)], [codar_line_lon(i); codar_line_lon(i+1)]);
end
mean_dist_codar = nanmean(dist_codar);

prof_per_codar = NaN(15,1);
diff_into_all = NaN(300,60);
diff_out_all = NaN(300,60);
binary_into_all = {};
binary_out_all = {};
bloom_binary_all = [];
threshold_all = [];
start_line_time = [];
bloom_binary_all_timestamp = [];
dens_diff_all = [];
bloom_binary_all_lon = [];
bloom_binary_all_lat = [];
rpd_all = [];
ftle_all = [];
mlparticle_all = [];

ad_poly_lon = [-64.02, -64.21, -64.3, -64.1, -64.01];
ad_poly_lat = [-64.82, -64.79, -64.86, -64.89, -64.82];
ge_poly_lon = [-63.85, -64.11, -63.95, -63.7, -63.8, -63.85];
ge_poly_lat = [-64.79, -64.89, -64.93, -64.86, -64.8, -64.79];

%% index for ADELIE transect and define phyto patches for each survey day
for i = 1:length(survey_days)
    day = survey_days(i,:);
    ind = find(ACRO.mtime >= datenum(day) & ACRO.mtime <= datenum(day)+1);
    
    lat = ACRO.lat(ind);
    lon = ACRO.lon(ind);
    dist = NaN(size(lon));
    mtime = ACRO.mtime(ind);
    dens_diff = ACRO.dens_diff(ind);
    chl = ACRO.mlparticle(ind); %%% edit here to change patch-determining variable
    rpd = ACRO.rpd_match(ind);
    ftle = ACRO.ftle_match(ind);
    mlparticle = ACRO.mlparticle(ind);
%     chl = max_part_all(ind);
    
    % index for just Adelie
    ind_ad = inpolygon(lon,lat,ad_poly_lon,ad_poly_lat);
    lat = lat(ind_ad);
    lon = lon(ind_ad);
    chl = chl(ind_ad);
    rpd = rpd(ind_ad);
    ftle = ftle(ind_ad);
    mtime = mtime(ind_ad);
    dens_diff = dens_diff(ind_ad);
    mlparticle = mlparticle(ind_ad);
    

    for j = 1:length(lat)-1
        [dist(j), AF, AR] = dist_lat_lon([lat(j); lat(j+1)], [lon(j); lon(j+1)]);
    end

    mean_dist_acro = nanmean(dist);
    prof_per_codar(i) = mean_dist_codar / mean_dist_acro;
    mean_dist_acro_all(i) = nanmean(dist);
%         bin_size = round(prof_per_codar(i));
    bin_size = 3;
    bin_plus1 = bin_size +1;
    bin_minus1 = bin_size -1;
    bin_minus2 = bin_size -2;
    bin_plus2 = bin_size +2;
    bin_minus3 = bin_size - 3;
    even_odd = rem(bin_size,2);
    
    % define threshold of patch ID with a sliding mean of bin size --> smooths ACRO data to resolution of HFR
    smoothd_chl =smoothdata(chl,'movmean',bin_size);
    med = nanmedian(chl);
    threshold = med*1.05;
    bloom = NaN(size(chl));

    if bin_size == 3
        for k = bin_minus1:length(chl)-(bin_minus1)
            bin = nanmean(chl(k-bin_minus2:1:k+bin_minus2));
            if bin >= threshold
                bloom(k) = 1;
            else
                bloom(k) = 0;
            end
        end
    elseif bin_size == 4
        for k = bin_size:length(chl)-(bin_minus3)
            bin = nanmean(chl(k-bin_minus2:1:k+bin_minus3));
            if bin >= threshold
                bloom(k) = 1;
            else
                bloom(k) = 0;
            end
        end
    else
        for k = bin_size:length(chl)-(bin_minus3)
            bin = nanmean(chl(k-bin_minus1:1:k));
            if bin >= threshold
                bloom(k) = 1;
            else
                bloom(k) = 0;
            end
        end
    end




    ind_bloom = find(bloom ==1);
    bloom_binary_all = [bloom_binary_all, bloom];
    bloom_binary_all_timestamp = [bloom_binary_all_timestamp, mtime];
    dens_diff_all = [dens_diff_all, dens_diff];
    bloom_binary_all_lon = [bloom_binary_all_lon, lon];
    bloom_binary_all_lat = [bloom_binary_all_lat, lat];
    rpd_all = [rpd_all, rpd];
    ftle_all = [ftle_all, ftle];
    mlparticle_all = [mlparticle_all, mlparticle];

    threshold_all = [threshold_all, threshold];
        

    
end


% save('patch_binary_smlparticle_max_linebyline.mat', 'bloom_binary_all');


%% index for ADELIE transect and define patch edges for each survey day
edges_all = [];
edges_threshold_all = [];

for i = 1:length(survey_days)
    day = survey_days(i,:);
    ind = find(ACRO.mtime >= datenum(day) & ACRO.mtime <= datenum(day)+1);
    
    lat = ACRO.lat(ind);
    lon = ACRO.lon(ind);
    dist = NaN(size(lon));
    mtime = ACRO.mtime(ind);
    chl = ACRO.mlparticle(ind); %%% edit here to change patch-determining variable
    
    % index for just Adelie
    ind_ad = inpolygon(lon,lat,ad_poly_lon,ad_poly_lat);
    lat = lat(ind_ad);
    lon = lon(ind_ad);
    chl = chl(ind_ad);
    mtime = mtime(ind_ad);
    

    for j = 1:length(lat)-1
        [dist(j), AF, AR] = dist_lat_lon([lat(j); lat(j+1)], [lon(j); lon(j+1)]);
    end

    mean_dist_acro = nanmean(dist);
    prof_per_codar(i) = mean_dist_codar / mean_dist_acro;
    mean_dist_acro_all(i) = nanmean(dist);
%         bin_size = round(prof_per_codar(i));
    bin_size = 3;
    bin_plus1 = bin_size +1;
    bin_minus1 = bin_size -1;
    bin_minus2 = bin_size -2;
    bin_plus2 = bin_size +2;
    bin_minus3 = bin_size - 3;
    even_odd = rem(bin_size,2);
    
    % define threshold of patch ID with a sliding mean of bin size --> smooths ACRO data to resolution of HFR
    smoothd_chl =smoothdata(chl,'movmean',bin_size);
    med = nanmedian(smoothd_chl);
    threshold = med*1.05;
    edges = NaN(size(chl));
    difference = NaN(size(chl));
    edges_threshold_diff = NaN(size(chl));
        
    %define threshold of change in chl
    diff_A = diff(smoothd_chl);
    diff_A = [NaN, diff_A];
    threshold_diff = nanmedian(abs(diff_A))*1.05;
        
        for k = (bin_size - 1):length(chl)-(bin_size - 1)
                j = k+1;
                difference(k) = abs(smoothd_chl(k)-smoothd_chl(j));
                if difference(k) >= threshold_diff
                    edges_threshold_diff(k) = 1;
                else
                    edges_threshold_diff(k) = 0;
                end
                if smoothd_chl(k) >= threshold && smoothd_chl(j) < threshold
                    edges(k) = 1;
                elseif smoothd_chl(k) >= threshold && smoothd_chl(j) >=threshold
                    edges(k) = 0;
                elseif smoothd_chl(k) < threshold && smoothd_chl(j) >=threshold
                    edges(k) = 2;
                elseif smoothd_chl(k) < threshold && smoothd_chl(j) < threshold
                    edges(k) = 0;
                else
                    edges(k) = NaN;
                end
        end
   

    edges_all = [edges_all, edges];
    edges_threshold_all = [edges_threshold_all, edges_threshold_diff];

    threshold_all = [threshold_all, threshold];
        

    
end


%% create patchID struct

patchID_Adelie.variable = 'integrated mixed layer particle backscatter';
patchID_Adelie.creationDate = '2-Nov-2022';
patchID_Adelie.createdWith = 'defining_bloom_final_reprocessed_adelie_transect.m';
patchID_Adelie.patch = bloom_binary_all;
patchID_Adelie.timestamp = bloom_binary_all_timestamp;
patchID_Adelie.dens_diff = dens_diff_all;
patchID_Adelie.lon = bloom_binary_all_lon;
patchID_Adelie.lat = bloom_binary_all_lat;
patchID_Adelie.ftle = ftle_all;
patchID_Adelie.rpd = rpd_all;
patchID_Adelie.edges = edges_all;
patchID_Adelie.edges_threshold = edges_threshold_all;
patchID_Adelie.mlparticle  = mlparticle_all;
