rng(0);
noiseVariancedBm = -96;
%center frequency of 2GHz
fc = 2e9;
lambda = 3e8/fc;
N = 25;
num_signal_snapshots = 200;
array_spacing = 0.5; %considering the antenna elements are spacedapart by lambda/2
half_angular_spread = 10; %20 degree angular spread
ula = phased.ULA('NumElements',N,'ElementSpacing',lambda/2);
sigma_sf = db2pow(8) %8 shadowing std of 8db 
aoa_sigma_offline = 2
kernel_function = "squaredexponential"

%Height of BS (in meters)
h_BS = 10;
%Height of UT (in meters)
h_UT = 1.5;
decorr = 13; %shadowing decorrelation distance

squareLength = 200; %each side of a square simulation area 
nbrOfSetups = 100; %this controls the number of AP-RP layouts i.e. number of setups with random UE and AP locations
num_tp_points = 1000; %consider 1000 random points for testing
L = 25  %total number of APs in simulation area
possible_RP_pos_per_row = [4,6,8,10,12,15];
%minimum distance between BSs and UEs
minDistanceUE2AP = 5.5; % (d2D > 5.5) => (d3D > 10m) for h_UT = 1m and h_BS = 10m; PL distance equation defined in Fraunhofer Region 
minDistanceAP2AP = 25; %neighbouring APs are spaced far apart to ensure the shadowing between two APs is uncorrelated

%Total uplink transmit power per UE (mW)
p = 100;

%comments correspond to a 25RP, 9AP setup
% Go through all RP positions
% For fingerprint database building, only one UE is placed at the RP location and RSS is measured.  

positioning_err_centr_rss = zeros(length(possible_RP_pos_per_row),nbrOfSetups,num_tp_points);
positioning_err_centr_aoa = zeros(length(possible_RP_pos_per_row),nbrOfSetups,num_tp_points);
positioning_err_centr_hybrid = zeros(length(possible_RP_pos_per_row),nbrOfSetups,num_tp_points);
positioning_err_distr_median = zeros(length(possible_RP_pos_per_row),nbrOfSetups,num_tp_points);
positioning_err_distr_mean = zeros(length(possible_RP_pos_per_row),nbrOfSetups,num_tp_points);
positioning_err_distr_zscore_tz_1 = zeros(length(possible_RP_pos_per_row),nbrOfSetups,num_tp_points);
positioning_err_distr_bayesian = zeros(length(possible_RP_pos_per_row),nbrOfSetups,num_tp_points);
positioning_err_distr_zscore_tz_1_bayesian = zeros(length(possible_RP_pos_per_row),nbrOfSetups,num_tp_points);

TP_positions_all_setup = cell(1, nbrOfSetups);
AP_positions_all_setup = cell(1, nbrOfSetups);

for setup_idx = 1:nbrOfSetups
    % Generate the setup and only store TP and AP positions
    [~, TP_positions_all_setup{setup_idx}, AP_positions_all_setup{setup_idx}] = ...
        cell_free_layout_setup_grid_RP(possible_RP_pos_per_row(end), ...
                                       squareLength, ...
                                       L, ...
                                       minDistanceUE2AP, ...
                                       minDistanceAP2AP, ...
                                       num_tp_points);
end

for RP_pos_idx = 1:length(possible_RP_pos_per_row)
    RP_positions_per_row = possible_RP_pos_per_row(RP_pos_idx)
    K = RP_positions_per_row^2;
    beta_fngprnt = zeros(K,L); %25x9
    nom_azi_angle_fngprnt = zeros(K,L);
    for setup_idx = 1:nbrOfSetups
        setup_idx
        [RP_positions,~,~] = cell_free_layout_setup_grid_RP(RP_positions_per_row,squareLength,L,minDistanceUE2AP,minDistanceAP2AP,num_tp_points);
        TP_positions = TP_positions_all_setup{setup_idx};
        AP_positions = AP_positions_all_setup{setup_idx};
        %functionPlotSetup(squareLength,RP_positions,AP_positions,TP_positions);
        %Prepare to store shadowing correlation matrix
        shadowCorrMatrix = sigma_sf^2*ones(K,K);
        shadowAPrealizations = zeros(K,L);
        for RP_idx = 1:K
            if RP_idx-1>0
                %Compute distances from the new prospective UE to all other UEs
                RPDistances = zeros(RP_idx-1,1);
                for i = 1:RP_idx-1
                    RPDistances(i) = abs(RP_positions(RP_idx) - RP_positions(i));
                end
                newcolumn = sigma_sf^2*2.^(-RPDistances/decorr);
                term1 = newcolumn'/shadowCorrMatrix(1:RP_idx-1,1:RP_idx-1);
                meanvalues = term1*shadowAPrealizations(1:RP_idx-1,:);
                stdvalue = sqrt(sigma_sf^2 - term1*newcolumn);
            else
                %Add the UE and begin to store shadow fading correlation values
                meanvalues = 0;
                stdvalue = sigma_sf;
                newcolumn = [];
            end
            %Generate the shadow fading realizations
            rss_shadowing_offline = meanvalues + stdvalue*randn(1,L);
            shadowCorrMatrix(1:RP_idx-1,RP_idx) = newcolumn;
            shadowCorrMatrix(RP_idx,1:RP_idx-1) = newcolumn';
            shadowAPrealizations(RP_idx,:) = rss_shadowing_offline;
            for AP_idx = 1:L
                aoa_err_offline = aoa_sigma_offline*randn;
                %disp(['Running RP' num2str(RP_idx) ' and AP ' num2str(AP_idx)]);
                d_2D = abs(RP_positions(RP_idx) - AP_positions(AP_idx));
                d_3D = sqrt((h_BS-h_UT)^2 + d_2D^2);     
                PL = 35.3*log10(d_3D) + 22.4 + 21.3*log10(fc/1e9);
                nom_azi_angle_fngprnt(RP_idx,AP_idx) = rad2deg(angle(RP_positions(RP_idx)- AP_positions(AP_idx))) + aoa_err_offline; %azimuth is angle between UE and BS on x axis
                beta_fngprnt(RP_idx,AP_idx) = -PL + rss_shadowing_offline(AP_idx);
            end %for AP_idx = 1:L 
        end %for RP_idx = 1:numel(RP_positions)
    
    finite_sample_effect_offline = gamrnd(num_signal_snapshots,1/num_signal_snapshots,K,L);
    RSS_fngprnt_mW = N*p*db2pow(beta_fngprnt) + N*db2pow(noiseVariancedBm);
    RSS_fngprnt_dB = 10*log10((RSS_fngprnt_mW.*finite_sample_effect_offline)/100); %25x9
    
    beta_fngprnt_tp = zeros(num_tp_points,L); %10000x9 matrix
    nom_azi_angle_TP_est = zeros(num_tp_points,L); %estimated azimuth angle by BS/CPU
    for TP_idx = 1:num_tp_points
        TPDistances = abs(RP_positions(:) - TP_positions(TP_idx)); %distance of the TP from each of the 25 RPs; 25x1 vector
        newcolumn = sigma_sf^2*2.^(-TPDistances/decorr);
        term1 = newcolumn'/shadowCorrMatrix;
        tp_meanvalues = term1*shadowAPrealizations;
        tp_stdvalue = sqrt(sigma_sf^2 - term1*newcolumn);
        rss_shadowing_online = tp_meanvalues + tp_stdvalue*randn(1,L);
        for AP_idx = 1:L
            d_2D_TP = abs(TP_positions(TP_idx)- AP_positions(AP_idx));
            d_3D_TP = sqrt((h_BS-h_UT)^2 + d_2D_TP^2);   
            theta = rad2deg(angle(TP_positions(TP_idx)- AP_positions(AP_idx)));
            broadside_angle = az2broadside(theta); %calculating broadside; This restricts angle to [-90 to 90]
            signal_pow_per_ant = p*db2pow(beta_fngprnt_tp(TP_idx,AP_idx))*(1e-3);
            noise_pow = db2pow(noiseVariancedBm)*(1e-3);
            [~,~,R] = custom_sensor_sig(getElementPosition(ula)/lambda,num_signal_snapshots,broadside_angle,half_angular_spread,array_spacing,noise_pow,signal_pow_per_ant);
            music_est = musicdoa(R,1,'ScanAngles',[-90:.1:90]);%R = covmat, 1 = 1 angle to be estimated
            if(theta > 90) %180degree ambiguity is resolved by considering for example, a second ULA at an angle to the main ULA
                nom_azi_angle_TP_est(TP_idx,AP_idx) = 180 - music_est;
            elseif(theta<-90)
                nom_azi_angle_TP_est(TP_idx,AP_idx) = -180 - music_est;
            else
                nom_azi_angle_TP_est(TP_idx,AP_idx) = music_est;
            end
            PL   = 35.3*log10(d_3D_TP) + 22.4 + 21.3*log10(fc/1e9);
            beta_fngprnt_tp(TP_idx,AP_idx) = -PL + rss_shadowing_online(AP_idx);
        end %for AP_idx = 1:L
    end %TP_idx = 1:num_tp_points
    
    finite_sample_effect_online = gamrnd(num_signal_snapshots,1/num_signal_snapshots,num_tp_points,L);
    RSS_tp_mW = N*p*db2pow(beta_fngprnt_tp) + N*db2pow(noiseVariancedBm); %RSS_tp_mW --> 16x9 vector
    RSS_tp_dB = 10*log10((RSS_tp_mW.*finite_sample_effect_online)/100); %16x9, with TPs in the centre of the 16 subregions
    
    %==========================================================================
    x_coordinates = real(RP_positions(:));
    y_coordinates = imag(RP_positions(:));
    
    RSS_var_names = strcat("RSS",num2str([1:L]'))';
    AoA_var_names = strcat("AoA",num2str([1:L]'))';

    var_names_x = [RSS_var_names AoA_var_names "x_coordinates"];
    var_names_y = [RSS_var_names AoA_var_names "y_coordinates"];
    
    centr_rss_gpr_X = array2table([RSS_fngprnt_dB,x_coordinates],"VariableNames",[RSS_var_names "x_coordinates"]);
    centr_rss_gpr_x_model = fitrgp(centr_rss_gpr_X,"x_coordinates",'kernelfunction',kernel_function);
    centr_rss_gpr_Y = array2table([RSS_fngprnt_dB,y_coordinates],"VariableNames",[RSS_var_names "y_coordinates"]);
    centr_rss_gpr_y_model = fitrgp(centr_rss_gpr_Y,"y_coordinates",'kernelfunction',kernel_function);

    centr_aoa_gpr_X = array2table([nom_azi_angle_fngprnt,x_coordinates],"VariableNames",[AoA_var_names "x_coordinates"]);
    centr_aoa_gpr_x_model = fitrgp(centr_aoa_gpr_X,"x_coordinates",'kernelfunction',kernel_function);
    centr_aoa_gpr_Y = array2table([nom_azi_angle_fngprnt,y_coordinates],"VariableNames",[AoA_var_names "y_coordinates"]);
    centr_aoa_gpr_y_model = fitrgp(centr_aoa_gpr_Y,"y_coordinates",'kernelfunction',kernel_function);
    
    centr_hybrid_gpr_X = array2table([RSS_fngprnt_dB,nom_azi_angle_fngprnt,x_coordinates],"VariableNames",var_names_x);
    centr_hybrid_gpr_x_model = fitrgp(centr_hybrid_gpr_X,"x_coordinates",'kernelfunction',kernel_function);
    centr_hybrid_gpr_Y = array2table([RSS_fngprnt_dB,nom_azi_angle_fngprnt,y_coordinates],"VariableNames",var_names_y);
    centr_hybrid_gpr_y_model = fitrgp(centr_hybrid_gpr_Y,"y_coordinates",'kernelfunction',kernel_function);

    for AP_idx = 1:L
        per_AP_RSS_fngprnt = RSS_fngprnt_dB(:,AP_idx);
        per_AP_AoA_fngprnt = nom_azi_angle_fngprnt(:,AP_idx);
        dyn_var_name_x = strcat('per_ap_gprx',num2str(AP_idx));
        dyn_var_name_y = strcat('per_ap_gpry',num2str(AP_idx));
        distr_gpr_X = table(per_AP_RSS_fngprnt,per_AP_AoA_fngprnt,x_coordinates);
        variable.(dyn_var_name_x) = fitrgp(distr_gpr_X,"x_coordinates",'kernelfunction',kernel_function);
        distr_gpr_Y = table(per_AP_RSS_fngprnt,per_AP_AoA_fngprnt,y_coordinates);
        variable.(dyn_var_name_y) = fitrgp(distr_gpr_Y,"y_coordinates",'kernelfunction',kernel_function);
    end

    %==========================================================================
    for TP_idx = 1:num_tp_points 
        per_tp_RSS = RSS_tp_dB(TP_idx,:); %per_tp_RSS --> 9x1 vector of the TP
        per_tp_est_angles = nom_azi_angle_TP_est(TP_idx,:);
        TRUE_TP_COORDS = [real(TP_positions(TP_idx)) imag(TP_positions(TP_idx))];
    
        [x_coordinate_tp] = predict(centr_rss_gpr_x_model, per_tp_RSS);
        [y_coordinate_tp] = predict(centr_rss_gpr_y_model, per_tp_RSS);
        CENTR_RSS_GPR_PRED = [x_coordinate_tp y_coordinate_tp];
        dist_coords = [CENTR_RSS_GPR_PRED;TRUE_TP_COORDS];
        positioning_err_centr_rss(RP_pos_idx,setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
    
        [x_coordinate_tp] = predict(centr_aoa_gpr_x_model, per_tp_est_angles);
        [y_coordinate_tp] = predict(centr_aoa_gpr_y_model, per_tp_est_angles);
        CENTR_AOA_GPR_PRED = [x_coordinate_tp y_coordinate_tp];
        dist_coords = [CENTR_AOA_GPR_PRED;TRUE_TP_COORDS];
        positioning_err_centr_aoa(RP_pos_idx,setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
    
        [x_coordinate_tp] = predict(centr_hybrid_gpr_x_model, [per_tp_RSS per_tp_est_angles]);
        [y_coordinate_tp] = predict(centr_hybrid_gpr_y_model, [per_tp_RSS per_tp_est_angles]);
        CENTR_HYBRID_GPR_PRED = [x_coordinate_tp y_coordinate_tp];
        dist_coords = [CENTR_HYBRID_GPR_PRED;TRUE_TP_COORDS];
        positioning_err_centr_hybrid(RP_pos_idx,setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
    
        DISTR_PRED_TP_X_COORDS = zeros(1,L); 
        DISTR_PRED_TP_Y_COORDS = zeros(1,L); 
        DISTR_PRED_TP_X_COORDS_SD = zeros(1,L); 
        DISTR_PRED_TP_Y_COORDS_SD = zeros(1,L);
        for AP_idx = 1:L
            AP_idx;
            dyn_var_name_x = strcat('per_ap_gprx',num2str(AP_idx));
            dyn_var_name_y = strcat('per_ap_gpry',num2str(AP_idx));
            [x_coordinate_tp,x_coordinate_sd] = predict(variable.(dyn_var_name_x), [per_tp_RSS(AP_idx) per_tp_est_angles(AP_idx)]);
            [y_coordinate_tp,y_coordinate_sd] = predict(variable.(dyn_var_name_y), [per_tp_RSS(AP_idx) per_tp_est_angles(AP_idx)]);
            DISTR_PRED_TP_X_COORDS(AP_idx) = x_coordinate_tp;
            DISTR_PRED_TP_Y_COORDS(AP_idx) = y_coordinate_tp; 
            DISTR_PRED_TP_X_COORDS_SD(AP_idx) = x_coordinate_sd;
            DISTR_PRED_TP_Y_COORDS_SD(AP_idx) = y_coordinate_sd; 
        end

        x_coordinate_median_tp = median(DISTR_PRED_TP_X_COORDS);
        y_coordinate_median_tp = median(DISTR_PRED_TP_Y_COORDS);
        DISTR_GPR_MEDIAN_PRED = [x_coordinate_median_tp y_coordinate_median_tp];
        dist_coords = [DISTR_GPR_MEDIAN_PRED;TRUE_TP_COORDS];
        positioning_err_distr_median(RP_pos_idx,setup_idx,TP_idx) = pdist(dist_coords,'euclidean');

        x_coordinate_mean_tp = mean(DISTR_PRED_TP_X_COORDS);
        y_coordinate_mean_tp = mean(DISTR_PRED_TP_Y_COORDS);
        DISTR_GPR_MEAN_PRED = [x_coordinate_mean_tp y_coordinate_mean_tp];
        dist_coords = [DISTR_GPR_MEAN_PRED;TRUE_TP_COORDS];
        positioning_err_distr_mean(RP_pos_idx,setup_idx,TP_idx) = pdist(dist_coords,'euclidean');

        x_coordinate_zscores = zscore(DISTR_PRED_TP_X_COORDS);
        y_coordinate_zscores = zscore(DISTR_PRED_TP_Y_COORDS);

        x_coordinate_zscore_mean_tp = mean(DISTR_PRED_TP_X_COORDS(abs(x_coordinate_zscores)<=1));
        y_coordinate_zscore_mean_tp = mean(DISTR_PRED_TP_Y_COORDS(abs(y_coordinate_zscores)<=1));
        DISTR_ZSCORE_GPR_PRED = [x_coordinate_zscore_mean_tp y_coordinate_zscore_mean_tp];
        dist_coords = [DISTR_ZSCORE_GPR_PRED;TRUE_TP_COORDS];
        positioning_err_distr_zscore_tz_1(RP_pos_idx,setup_idx,TP_idx) = pdist(dist_coords,'euclidean');

        x_coordinate_bayesian_var = 1 / sum(1 ./ (DISTR_PRED_TP_X_COORDS_SD .^ 2));
        x_coordinate_bayesian_sd = sqrt(x_coordinate_bayesian_var);
        x_coordinate_bayesian_tp = x_coordinate_bayesian_var * sum(DISTR_PRED_TP_X_COORDS ./ (DISTR_PRED_TP_X_COORDS_SD .^ 2));

        y_coordinate_bayesian_var = 1 / sum(1 ./ (DISTR_PRED_TP_Y_COORDS_SD .^ 2));
        y_coordinate_bayesian_sd = sqrt(y_coordinate_bayesian_var);
        y_coordinate_bayesian_tp = y_coordinate_bayesian_var * sum(DISTR_PRED_TP_Y_COORDS ./ (DISTR_PRED_TP_Y_COORDS_SD .^ 2));

        DISTR_BAY_GPR_PRED = [x_coordinate_bayesian_tp y_coordinate_bayesian_tp];
        dist_coords = [DISTR_BAY_GPR_PRED;TRUE_TP_COORDS];
        positioning_err_distr_bayesian(RP_pos_idx,setup_idx,TP_idx) = pdist(dist_coords,'euclidean');

        x_coordinate_zscores = zscore(DISTR_PRED_TP_X_COORDS);
        y_coordinate_zscores = zscore(DISTR_PRED_TP_Y_COORDS);

        % Find the indices where the absolute z-score is <= 1
        x_valid_indices = abs(x_coordinate_zscores) <= 1;
        y_valid_indices = abs(y_coordinate_zscores) <= 1;

        % Extract the corresponding mean values using logical indexing
        corresponding_means_x = DISTR_PRED_TP_X_COORDS(x_valid_indices);
        corresponding_means_y = DISTR_PRED_TP_Y_COORDS(y_valid_indices);

        % Extract the corresponding SD values using logical indexing
        corresponding_sds_x = DISTR_PRED_TP_X_COORDS_SD(x_valid_indices);
        corresponding_sds_y = DISTR_PRED_TP_Y_COORDS_SD(y_valid_indices);

        x_coordinate_zscore_bayesian_var = 1 / sum(1 ./ (corresponding_sds_x .^ 2));
        x_coordinate_zscore_bayesian_sd = sqrt(x_coordinate_zscore_bayesian_var);
        x_coordinate_zscore_bayesian_tp = x_coordinate_zscore_bayesian_var * sum(corresponding_means_x ./ (corresponding_sds_x .^ 2));

        y_coordinate_zscore_bayesian_var = 1 / sum(1 ./ (corresponding_sds_y .^ 2));
        y_coordinate_zscore_bayesian_sd = sqrt(y_coordinate_zscore_bayesian_var);
        y_coordinate_zscore_bayesian_tp = y_coordinate_zscore_bayesian_var * sum(corresponding_means_y ./ (corresponding_sds_y .^ 2));

        
        DISTR_BAY_ZSCORE_GPR_PRED = [x_coordinate_zscore_bayesian_tp y_coordinate_zscore_bayesian_tp];
        dist_coords = [DISTR_BAY_ZSCORE_GPR_PRED;TRUE_TP_COORDS];
        positioning_err_distr_zscore_tz_1_bayesian(RP_pos_idx,setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
    
    end
    fprintf("positioning_err_centr_rss = %d\n",mean(positioning_err_centr_rss(RP_pos_idx,setup_idx,:)));
    fprintf("positioning_err_centr_aoa = %d\n",mean(positioning_err_centr_aoa(RP_pos_idx,setup_idx,:)));
    fprintf("positioning_err_centr_hybrid = %d\n",mean(positioning_err_centr_rss(RP_pos_idx,setup_idx,:)));
    fprintf("positioning_err_distr_median = %d\n",mean(positioning_err_distr_median(RP_pos_idx,setup_idx,:)));
    fprintf("positioning_err_distr_mean = %d\n",mean(positioning_err_distr_mean(RP_pos_idx,setup_idx,:)));
    fprintf("positioning_err_distr_zscore_tz_1 = %d\n",mean(positioning_err_distr_zscore_tz_1(RP_pos_idx,setup_idx,:)));
    fprintf("positioning_err_distr_bayesian = %d\n",mean(positioning_err_distr_bayesian(RP_pos_idx,setup_idx,:)));
    fprintf("positioning_err_distr_zscore_tz_1_bayesian = %d\n",mean(positioning_err_distr_zscore_tz_1_bayesian(RP_pos_idx,setup_idx,:)));
    end %for n = 1:nbrOfSetups
end %RP_pos_idx = 1:length(possible_RP_pos_per_row)