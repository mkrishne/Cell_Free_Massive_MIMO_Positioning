rng(8);
diary Fig3b_localization_across_N_diary;

noiseVariancedBm = -96;
%center frequency of 2GHz
format long
fc = 2e9;
lambda = 3e8/fc;
% N = Number of antennas per AP
N = [3,4,5,6,8,10,12,15,17,19,21,25,30,35,40,45,50];
N_max = length(N); %sweep from 2 antenna per AP to 100 antenna per AP
array_spacing = 0.5;
half_angular_spread = 10; %20 degree angular spread
num_signal_snapshots = 200;
sigma_sf = db2pow(8) %8 shadowing std of 8db 
aoa_sigma_offline = 2
kernel_function = "squaredexponential"

%Height of BS (in meters)
h_BS = 10;
%Height of UT (in meters)
h_UT = 1.5;
decorr = 13; %shadowing decorrelation distance
SIGMA_BAR = 2;

squareLength = 200; %each side of a square simulation area 
nbrOfSetups = 100; %this controls the number of AP-RP layouts i.e. number of setups with random UE and AP locations
num_tp_points = 1000; %consider 1000 random points in the simulation area for testing
L = 25 %total number of APs in simulation area
RP_positions_per_row = 15 %creates a grid of (RP_positions_per_row x RP_positions_per_row) RP positions
%RP_positions_per_row = 8 %for FIg3a, K=64, set RP_positions_per_row=8
K = RP_positions_per_row^2;
% Number of nearest neighbors
KNN_K = 4;

%minimum distance between BSs and UEs
minDistanceUE2AP = 5.5; % (d2D > 5.5) => (d3D > 10m) for h_UT = 1m and h_BS = 10m; PL distance equation defined in Fraunhofer Region 
minDistanceAP2AP = 25; %neighbouring APs are spaced far apart to ensure the shadowing between two APs is uncorrelated

%Total uplink transmit power per UE (mW)
p = 100;

%comments correspond to a 25RP, 9AP setup
% Go through all RP positions
% For fingerprint database building, only one UE is placed at the RP location and RSS is measured.  

avg_positioning_err_centr_rss = zeros(1,N_max);
avg_positioning_err_centr_aoa = zeros(1,N_max);
avg_positioning_err_centr_hybrid = zeros(1,N_max);

avg_positioning_err_distr_median = zeros(1,N_max);
avg_positioning_err_distr_mean = zeros(1,N_max);
avg_positioning_err_distr_zscore_tz_0p6 = zeros(1,N_max);
avg_positioning_err_distr_zscore_tz_0p8 = zeros(1,N_max);
avg_positioning_err_distr_zscore_tz_1 = zeros(1,N_max);
avg_positioning_err_distr_zscore_tz_1p2 = zeros(1,N_max);
avg_positioning_err_distr_zscore_tz_1p4 = zeros(1,N_max);

avg_positioning_err_distr_bayesian = zeros(1,N_max);
avg_positioning_err_distr_zscore_tz_0p6_bayesian = zeros(1,N_max);
avg_positioning_err_distr_zscore_tz_0p8_bayesian = zeros(1,N_max);
avg_positioning_err_distr_zscore_tz_1_bayesian = zeros(1,N_max);
avg_positioning_err_distr_zscore_tz_1p2_bayesian = zeros(1,N_max);
avg_positioning_err_distr_zscore_tz_1p4_bayesian = zeros(1,N_max);

avg_positioning_err_centr_hybrid_crlb = zeros(1,N_max);
avg_positioning_err_centr_aoa_crlb = zeros(1,N_max);
avg_positioning_err_distr_median_crlb = zeros(1,N_max);
avg_positioning_err_distr_mean_crlb = zeros(1,N_max);
avg_positioning_err_distr_bayesian_crlb = zeros(1,N_max);
avg_positioning_err_distr_zscore_tz_1_crlb = zeros(1,N_max);
avg_positioning_err_distr_zscore_tz_1_bayesian_crlb = zeros(1,N_max);
%-----------------------------------------------------------
avg_err_ellipse_area_centr_rss = zeros(1,N_max);
avg_err_ellipse_area_centr_aoa = zeros(1,N_max);
avg_err_ellipse_area_centr_hybrid = zeros(1,N_max);

avg_err_ellipse_area_distr_median = zeros(1,N_max);
avg_err_ellipse_area_distr_mean = zeros(1,N_max);
avg_err_ellipse_area_distr_bayesian = zeros(1,N_max);
avg_err_ellipse_area_distr_zscore_tz_0p6 = zeros(1,N_max);
avg_err_ellipse_area_distr_zscore_tz_0p8 = zeros(1,N_max);
avg_err_ellipse_area_distr_zscore_tz_1 = zeros(1,N_max);
avg_err_ellipse_area_distr_zscore_tz_1p2 = zeros(1,N_max);
avg_err_ellipse_area_distr_zscore_tz_1p4 = zeros(1,N_max);
avg_err_ellipse_area_distr_zscore_tz_0p6_bayesian = zeros(1,N_max);
avg_err_ellipse_area_distr_zscore_tz_0p8_bayesian = zeros(1,N_max);
avg_err_ellipse_area_distr_zscore_tz_1_bayesian = zeros(1,N_max);
avg_err_ellipse_area_distr_zscore_tz_1p2_bayesian = zeros(1,N_max);
avg_err_ellipse_area_distr_zscore_tz_1p4_bayesian = zeros(1,N_max);

avg_err_ellipse_area_centr_hybrid_crlb = zeros(1,N_max);
avg_err_ellipse_area_centr_aoa_crlb = zeros(1,N_max);
avg_err_ellipse_area_distr_median_crlb = zeros(1,N_max);
avg_err_ellipse_area_distr_mean_crlb = zeros(1,N_max);
avg_err_ellipse_area_distr_bayesian_crlb = zeros(1,N_max);
avg_err_ellipse_area_distr_zscore_tz_1_crlb = zeros(1,N_max);
avg_err_ellipse_area_distr_zscore_tz_1_bayesian_crlb = zeros(1,N_max);
%-----------------------------------------------------------
avg_positioning_err_centr_rss_LR = zeros(1,N_max);
avg_positioning_err_centr_aoa_LR = zeros(1,N_max);
avg_positioning_err_centr_hybrid_LR = zeros(1,N_max);

avg_positioning_err_centr_rss_LR_dB = zeros(1,N_max);
avg_positioning_err_centr_hybrid_LR_dB = zeros(1,N_max);
avg_positioning_err_distr_median_LR_dB = zeros(1,N_max);

avg_positioning_err_centr_aoa_LR_crlb = zeros(1,N_max);
avg_positioning_err_centr_hybrid_LR_crlb = zeros(1,N_max);
avg_positioning_err_centr_hybrid_LR_dB_crlb = zeros(1,N_max);
avg_positioning_err_distr_median_LR_dB_crlb = zeros(1,N_max);

avg_positioning_err_centr_rss_WKNN = zeros(1,N_max);
avg_positioning_err_centr_aoa_WKNN = zeros(1,N_max);
avg_positioning_err_centr_hybrid_WKNN = zeros(1,N_max);

avg_positioning_err_centr_rss_WKNN_dB = zeros(1,N_max);
avg_positioning_err_centr_hybrid_WKNN_dB = zeros(1,N_max);
avg_positioning_err_distr_median_WKNN_dB = zeros(1,N_max);

avg_positioning_err_centr_aoa_WKNN_crlb = zeros(1,N_max);
avg_positioning_err_centr_hybrid_WKNN_dB_crlb = zeros(1,N_max);
avg_positioning_err_distr_median_WKNN_dB_crlb = zeros(1,N_max);

avg_positioning_err_fcnn = zeros(1,N_max);
positioning_err_fcnn = zeros(nbrOfSetups,1);
%-----------------------------------------------------------
per_setup_positioning_err_centr_rss = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_centr_aoa = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_centr_hybrid = zeros(nbrOfSetups,num_tp_points);

per_setup_positioning_err_distr_median = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_mean = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_bayesian = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_zscore_tz_0p6 = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_zscore_tz_0p8 = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_zscore_tz_1 = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_zscore_tz_1p2 = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_zscore_tz_1p4 = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_zscore_tz_0p6_bayesian = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_zscore_tz_0p8_bayesian = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_zscore_tz_1_bayesian = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_zscore_tz_1p2_bayesian = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_zscore_tz_1p4_bayesian = zeros(nbrOfSetups,num_tp_points);

per_setup_positioning_err_centr_hybrid_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_centr_aoa_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_median_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_mean_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_bayesian_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_zscore_tz_1_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_zscore_tz_1_bayesian_crlb = zeros(nbrOfSetups,num_tp_points);

%-----------------------------------------------------------
per_setup_err_ellipse_area_centr_rss = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_centr_aoa = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_centr_hybrid = zeros(nbrOfSetups,num_tp_points);

per_setup_err_ellipse_area_distr_median = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_mean = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_bayesian = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_zscore_tz_0p6 = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_zscore_tz_0p8 = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_zscore_tz_1 = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_zscore_tz_1p2 = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_zscore_tz_1p4 = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_zscore_tz_0p6_bayesian = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_zscore_tz_0p8_bayesian = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_zscore_tz_1_bayesian = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_zscore_tz_1p2_bayesian = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_zscore_tz_1p4_bayesian = zeros(nbrOfSetups,num_tp_points);

per_setup_err_ellipse_area_centr_hybrid_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_centr_aoa_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_median_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_mean_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_bayesian_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_zscore_tz_1_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_err_ellipse_area_distr_zscore_tz_1_bayesian_crlb = zeros(nbrOfSetups,num_tp_points);

%-----------------------------------------------------------
per_setup_positioning_err_centr_rss_LR = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_centr_aoa_LR = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_centr_hybrid_LR = zeros(nbrOfSetups,num_tp_points);

per_setup_positioning_err_centr_rss_LR_dB = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_centr_hybrid_LR_dB = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_median_LR_dB = zeros(nbrOfSetups,num_tp_points);

per_setup_positioning_err_centr_aoa_LR_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_centr_hybrid_LR_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_centr_hybrid_LR_dB_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_median_LR_dB_crlb = zeros(nbrOfSetups,num_tp_points);

per_setup_positioning_err_centr_rss_WKNN = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_centr_aoa_WKNN = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_centr_hybrid_WKNN = zeros(nbrOfSetups,num_tp_points);

per_setup_positioning_err_centr_rss_WKNN_dB = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_centr_hybrid_WKNN_dB = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_median_WKNN_dB = zeros(nbrOfSetups,num_tp_points);

per_setup_positioning_err_centr_aoa_WKNN_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_centr_hybrid_WKNN_dB_crlb = zeros(nbrOfSetups,num_tp_points);
per_setup_positioning_err_distr_median_WKNN_dB_crlb = zeros(nbrOfSetups,num_tp_points);
%-----------------------------------------------------------

NUM_TP_INSIDE_ELLIPSE_CR = 0; %number of true TPs within the error ellipse for Centralized RSS
NUM_TP_INSIDE_ELLIPSE_CA = 0; %Centralized AOA
NUM_TP_INSIDE_ELLIPSE_CH = 0; %Centralized Hybrid
NUM_TP_INSIDE_ELLIPSE_DD = 0; %Distributed meDian
NUM_TP_INSIDE_ELLIPSE_DM = 0; %Distributed Mean
NUM_TP_INSIDE_ELLIPSE_DZ = 0; %Distributed Mean Zscore
NUM_TP_INSIDE_ELLIPSE_DBAY = 0; %Distributed Bayesian
NUM_TP_INSIDE_ELLIPSE_DZ_BAY = 0; %Distributed Bayesian Zscore

NUM_TP_INSIDE_ELLIPSE_CA_CRLB = 0; %Centralized AOA CRLB
NUM_TP_INSIDE_ELLIPSE_CH_CRLB = 0; %Centralized Hybrid CRLB
NUM_TP_INSIDE_ELLIPSE_DD_CRLB = 0; %Distributed meDian CRLB
NUM_TP_INSIDE_ELLIPSE_DM_CRLB = 0; %Distributed Mean CRLB
NUM_TP_INSIDE_ELLIPSE_DZ_CRLB = 0; %Distributed Zscore CRLB
NUM_TP_INSIDE_ELLIPSE_DBAY_CRLB = 0; %Distributed Bayesian CRLB
NUM_TP_INSIDE_ELLIPSE_DZ_BAY_CRLB = 0; %Distributed Bayesian Zscore CRLB

% Initialize cell arrays to store setups and shadowing data
RP_positions_all_setup = cell(1, nbrOfSetups);
TP_positions_all_setup = cell(1, nbrOfSetups);
AP_positions_all_setup = cell(1, nbrOfSetups);

for setup_idx = 1:nbrOfSetups
    % Generate the setup
    [RP_positions_all_setup{setup_idx}, ...
     TP_positions_all_setup{setup_idx}, ...
     AP_positions_all_setup{setup_idx}] = ...
     cell_free_layout_setup_random_RP(RP_positions_per_row, ...
                                      squareLength, ...
                                      L, ...
                                      minDistanceUE2AP, ...
                                      minDistanceAP2AP, ...
                                      num_tp_points);
end

shadowCorrMatrix_all_setup = cell(1, nbrOfSetups);
shadowAPrealizations_all_setup = cell(1, nbrOfSetups);
rss_shadowing_offline_all_setup = cell(1, nbrOfSetups);
for setup_idx = 1:nbrOfSetups
    RP_positions = RP_positions_all_setup{setup_idx};
    TP_positions = TP_positions_all_setup{setup_idx};
    AP_positions = AP_positions_all_setup{setup_idx};

    % Initialize shadowing correlation matrix and shadow fading realizations
    shadowCorrMatrix = sigma_sf^2 * ones(K, K);
    shadowAPrealizations = zeros(K, L);
    rss_shadowing_offline_list = cell(K, 1); % Store rss_shadowing_offline for each RP_idx

    for RP_idx = 1:K
        if RP_idx > 1
            %Compute distances from the new prospective UE to all other UEs
            RPDistances = zeros(RP_idx-1,1);
            for i = 1:RP_idx-1
                RPDistances(i) = abs(RP_positions(RP_idx) - RP_positions(i));
            end
            newcolumn = sigma_sf^2 * 2.^(-RPDistances / decorr);
            term1 = newcolumn' / shadowCorrMatrix(1:RP_idx-1, 1:RP_idx-1);
            meanvalues = term1 * shadowAPrealizations(1:RP_idx-1, :);
            stdvalue = sqrt(sigma_sf^2 - term1 * newcolumn);
        else
            % Add the UE and begin to store shadow fading correlation values
            meanvalues = 0;
            stdvalue = sigma_sf;
            newcolumn = [];
        end
        % Generate the shadow fading realizations
        rss_shadowing_offline = meanvalues + stdvalue * randn(1, L);
        shadowCorrMatrix(1:RP_idx-1, RP_idx) = newcolumn;
        shadowCorrMatrix(RP_idx, 1:RP_idx-1) = newcolumn';
        shadowAPrealizations(RP_idx, :) = rss_shadowing_offline;

        % Store rss_shadowing_offline for this RP_idx
        rss_shadowing_offline_list{RP_idx} = rss_shadowing_offline;
    end

    % Store shadowCorrMatrix, shadowAPrealizations, and rss_shadowing_offline for this setup
    shadowCorrMatrix_all_setup{setup_idx} = shadowCorrMatrix;
    shadowAPrealizations_all_setup{setup_idx} = shadowAPrealizations;
    rss_shadowing_offline_all_setup{setup_idx} = rss_shadowing_offline_list;
end

rss_shadowing_online_all_setup = cell(1, nbrOfSetups);
for setup_idx = 1:nbrOfSetups
    % Retrieve previously stored values for the current setup
    RP_positions = RP_positions_all_setup{setup_idx};
    TP_positions = TP_positions_all_setup{setup_idx};
    shadowCorrMatrix = shadowCorrMatrix_all_setup{setup_idx};
    shadowAPrealizations = shadowAPrealizations_all_setup{setup_idx};

    % Initialize storage for rss_shadowing_online for this setup
    rss_shadowing_online_list = cell(1, num_tp_points);

    % Loop over all test points (TP_idx)
    for TP_idx = 1:num_tp_points
        % Compute distances between reference points and the current test point
        TPDistances = abs(RP_positions(:) - TP_positions(TP_idx));
        newcolumn = sigma_sf^2 * 2.^(-TPDistances / decorr);
        term1 = newcolumn' / shadowCorrMatrix;
        tp_meanvalues = term1 * shadowAPrealizations;
        tp_stdvalue = sqrt(sigma_sf^2 - term1 * newcolumn);

        % Generate shadowing realizations for the current test point
        rss_shadowing_online = tp_meanvalues + tp_stdvalue * randn(1, L);

        % Store rss_shadowing_online for this test point
        rss_shadowing_online_list{TP_idx} = rss_shadowing_online;
    end

    % Store rss_shadowing_online for this setup
    rss_shadowing_online_all_setup{setup_idx} = rss_shadowing_online_list;
end

for n_idx = 1:N_max
    ant_number = N(n_idx)
    ula = phased.ULA('NumElements',N(n_idx),'ElementSpacing',lambda/2);
    puncture_idx = [];
    for setup_idx = 1:nbrOfSetups
        setup_idx
        beta_fngprnt = zeros(K,L); %25x9
        nom_azi_angle_fngprnt = zeros(K,L);
        psi_fngprnt = zeros(K,L); %25x9
        F = dftmtx(N(n_idx))/sqrt(N(n_idx));
        FH = F';  % Hermitian transpose of F

        RP_positions = RP_positions_all_setup{setup_idx};
        TP_positions = TP_positions_all_setup{setup_idx};
        AP_positions = AP_positions_all_setup{setup_idx};
        %functionPlotSetup(squareLength,RP_positions,AP_positions,TP_positions);
        %Prepare to store shadowing correlation matrix
        rss_shadowing_offline_list = rss_shadowing_offline_all_setup{setup_idx};
        rss_shadowing_online_list = rss_shadowing_online_all_setup{setup_idx};

        for RP_idx = 1:K
            rss_shadowing_offline = rss_shadowing_offline_list{RP_idx};
            for AP_idx = 1:L
                aoa_err_offline = aoa_sigma_offline*randn;
                %disp(['Running RP' num2str(RP_idx) ' and AP ' num2str(AP_idx)]);
                d_2D = abs(RP_positions(RP_idx) - AP_positions(AP_idx));
                d_3D = sqrt((h_BS-h_UT)^2 + d_2D^2);     
                PL = 35.3*log10(d_3D) + 22.4 + 21.3*log10(fc/1e9);
                beta_fngprnt(RP_idx,AP_idx) = -PL + rss_shadowing_offline(AP_idx);
                nom_azi_angle = rad2deg(angle(RP_positions(RP_idx)- AP_positions(AP_idx)));
                nom_azi_angle_fngprnt(RP_idx,AP_idx) = nom_azi_angle + aoa_err_offline; %azimuth is angle between UE and BS on x axis
                signal_pow_per_ant = p*db2pow(beta_fngprnt(RP_idx,AP_idx))*(1e-3); %power per antenna in watt
                noise_pow = db2pow(noiseVariancedBm)*(1e-3); %noise power in watt
                [~,~,C,~] = custom_sensor_sig(getElementPosition(ula)/lambda,num_signal_snapshots,nom_azi_angle,half_angular_spread,array_spacing,noise_pow,signal_pow_per_ant);
                FHCF = FH * C * F;  % Compute F^H C F
                psi_fngprnt(RP_idx,AP_idx) = real(trace(FHCF));  % Calculate the trace
            end %for AP_idx = 1:L
        end %for RP_idx = 1:numel(RP_positions)
    
    finite_sample_effect_offline = gamrnd(num_signal_snapshots,1/num_signal_snapshots,K,L);
    RSS_fngprnt_mW = N(n_idx)*p*db2pow(beta_fngprnt) + N(n_idx)*db2pow(noiseVariancedBm);
    RSS_fngprnt_dB = 10*log10((RSS_fngprnt_mW.*finite_sample_effect_offline)/100); %25x9
    psi_fngprnt_dB = 10*log10(psi_fngprnt*1e3/100); %25x9

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ONLINE STAGE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta_fngprnt_tp = zeros(num_tp_points,L); %10000x9 matrix
    psi_fngprnt_tp = zeros(num_tp_points,L); %10000x9 matrix
    nom_azi_angle_TP_music_est = zeros(num_tp_points,L); %estimated azimuth angle by BS/CPU using MUSIC algo
    nom_azi_angle_TP_crlb = zeros(num_tp_points,L);
    
    for TP_idx = 1:num_tp_points
        rss_shadowing_online = rss_shadowing_online_list{TP_idx};
        for AP_idx = 1:L
            d_2D_TP = abs(TP_positions(TP_idx)- AP_positions(AP_idx));
            d_3D_TP = sqrt((h_BS-h_UT)^2 + d_2D_TP^2);   
            PL   = 35.3*log10(d_3D_TP) + 22.4 + 21.3*log10(fc/1e9);
            beta_fngprnt_tp(TP_idx,AP_idx) = -PL + rss_shadowing_online(AP_idx);
            theta = rad2deg(angle(TP_positions(TP_idx)- AP_positions(AP_idx)));
            broadside_angle = az2broadside(theta); %calculating broadside; This restricts angle to [-90 to 90]
            %In the az2broadside calculation above, the array's normal axis is aligned with the x-axis to facilitate the use of MATLAB's convenient built-in functions.
            signal_pow_per_ant = p*db2pow(beta_fngprnt_tp(TP_idx,AP_idx))*(1e-3); %power per antenna in watt
            noise_pow = db2pow(noiseVariancedBm)*(1e-3); %noise power in watt
            [~,R_theoretical,R,steeringvec] = custom_sensor_sig(getElementPosition(ula)/lambda,num_signal_snapshots,broadside_angle,half_angular_spread,array_spacing,noise_pow,signal_pow_per_ant);
            music_est = musicdoa(R,1,'ScanAngles',[-90:.1:90]);%R = covmat, 1 = 1 angle to be estimated
            if(theta > 90) %180degree ambiguity is resolved by considering for example, a second ULA at an angle to the main ULA or analysis of the estimated angles of multiple APs
                nom_azi_angle_TP_music_est(TP_idx,AP_idx) = 180 - music_est;
            elseif(theta<-90)
                nom_azi_angle_TP_music_est(TP_idx,AP_idx) = -180 - music_est;
            else
                nom_azi_angle_TP_music_est(TP_idx,AP_idx) = music_est;
            end
            
            crlb = find_crlb(num_signal_snapshots,array_spacing,broadside_angle,half_angular_spread,steeringvec,signal_pow_per_ant,R);
            aoa_err_crlb = rad2deg(sqrt(crlb))*randn;
            nom_azi_angle_TP_crlb(TP_idx,AP_idx) = theta + aoa_err_crlb; %azimuth is angle between UE and BS on x axis
            
            FHCtpF = FH * R * F;  % Compute F^H C F
            psi_fngprnt_tp(TP_idx,AP_idx) = real(trace(FHCtpF));  % Calculate the trace
         end %for AP_idx = 1:L
    end %TP_idx = 1:num_tp_points
    
    finite_sample_effect_online = gamrnd(num_signal_snapshots,1/num_signal_snapshots,num_tp_points,L);
    RSS_tp_mW = N(n_idx)*p*db2pow(beta_fngprnt_tp) + N(n_idx)*db2pow(noiseVariancedBm); %RSS_tp_mW --> 16x9 vector
    RSS_tp_dB = 10*log10((RSS_tp_mW.*finite_sample_effect_online)/100); %16x9, with TPs in the centre of the 16 subregions
    psi_fngprnt_tp_dB = 10*log10(psi_fngprnt_tp*1e3/100); %25x9

    %============================== ML MODEL TRAINING ============================================
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

    y_train_LR = [real(RP_positions(:)), imag(RP_positions(:))]; % RP positions (real and imaginary parts as separate columns)
    X_train_LR = [ones(K,1) psi_fngprnt]; % Feature matrix for RPs
    X_train_LR_AOA = [ones(K,1) nom_azi_angle_fngprnt]; % Feature matrix for RPs
    X_train_LR_hybrid = [ones(K,1) psi_fngprnt nom_azi_angle_fngprnt]; % Feature matrix for RPs

    X_train_LR_dB = [ones(K,1) psi_fngprnt_dB]; % Feature matrix for RPs
    X_train_LR_hybrid_dB = [ones(K,1) psi_fngprnt_dB nom_azi_angle_fngprnt]; % Feature matrix for RPs

    % Train model for x and y coordinates separately
    w_x = regress(y_train_LR(:,1), X_train_LR);
    w_y = regress(y_train_LR(:,2), X_train_LR);
    w_x_aoa = regress(y_train_LR(:,1), X_train_LR_AOA);
    w_y_aoa = regress(y_train_LR(:,2), X_train_LR_AOA);
    w_x_hybrid = regress(y_train_LR(:,1), X_train_LR_hybrid);
    w_y_hybrid = regress(y_train_LR(:,2), X_train_LR_hybrid);
   
    w_x_dB = regress(y_train_LR(:,1), X_train_LR_dB);
    w_y_dB = regress(y_train_LR(:,2), X_train_LR_dB);
    w_x_hybrid_dB = regress(y_train_LR(:,1), X_train_LR_hybrid_dB);
    w_y_hybrid_dB = regress(y_train_LR(:,2), X_train_LR_hybrid_dB);

    y_train_KNN = RP_positions;
    X_train_KNN = psi_fngprnt; % Feature matrix for RPs        
    X_train_KNN_AOA = nom_azi_angle_fngprnt; % Feature matrix for RPs
    X_train_KNN_hybrid = [psi_fngprnt nom_azi_angle_fngprnt]; % Feature matrix for RPs

    X_train_KNN_dB = psi_fngprnt_dB; % Feature matrix for RPs 
    X_train_KNN_hybrid_dB = [psi_fngprnt_dB nom_azi_angle_fngprnt]; % Feature matrix for RPs

    X_train_KNN_AOA_crlb = nom_azi_angle_fngprnt; % Feature matrix for RPs
    X_train_KNN_hybrid_dB_crlb = [psi_fngprnt_dB nom_azi_angle_fngprnt]; % Feature matrix for RPs

    for AP_idx = 1:L
        per_AP_RSS_fngprnt = psi_fngprnt_dB(:,AP_idx);
        per_AP_AoA_fngprnt = nom_azi_angle_fngprnt(:,AP_idx);
        X_train_LR_hybrid_per_ap = [ones(K,1) per_AP_RSS_fngprnt per_AP_AoA_fngprnt]; % per_ap means distributed
        dyn_var_name_x = strcat('w_per_ap_x',num2str(AP_idx));
        dyn_var_name_y = strcat('w_per_ap_y',num2str(AP_idx));
        variable.(dyn_var_name_x) = regress(y_train_LR(:,1), X_train_LR_hybrid_per_ap);
        variable.(dyn_var_name_y) = regress(y_train_LR(:,2), X_train_LR_hybrid_per_ap);

        dyn_var_name_knn = strcat('X_train_KNN_hybrid_per_ap',num2str(AP_idx));
        variable.(dyn_var_name_knn) = [per_AP_RSS_fngprnt per_AP_AoA_fngprnt];
    end

    X_NN = nom_azi_angle_fngprnt';
    Y_NN = [x_coordinates, y_coordinates]';
    hiddenLayerSizes = [128, 64, 32, 32, 16];
    fcnn_model = feedforwardnet(hiddenLayerSizes, 'trainscg');

    % Configure training parameters
    fcnn_model.trainParam.epochs = 500;
    fcnn_model.trainParam.goal = 1e-6;
    fcnn_model.trainParam.min_grad = 1e-6;
    fcnn_model.trainParam.showWindow = false;

    fcnn_model.divideParam.trainRatio = 1.0;
    fcnn_model.divideParam.valRatio   = 0.0;
    fcnn_model.divideParam.testRatio  = 0.0;
    
    % Train the network
    fcnn_model = train(fcnn_model, X_NN, Y_NN, 'useGPU', 'yes');
    %==========================================================================
    for TP_idx = 1:num_tp_points 
        per_tp_RSS = RSS_tp_dB(TP_idx,:); %per_tp_RSS --> 9x1 vector of the TP
        per_tp_est_angles = nom_azi_angle_TP_music_est(TP_idx,:);
        per_tp_est_angles_crlb = nom_azi_angle_TP_crlb(TP_idx,:);
        TRUE_TP_COORDS = [real(TP_positions(TP_idx)) imag(TP_positions(TP_idx))];
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%% CENTRALISED %%%%%%%%%%%%%%%%%%%%%%%%%%%
        [x_coordinate_tp,x_coordinate_sd] = predict(centr_rss_gpr_x_model, per_tp_RSS);
        [y_coordinate_tp,y_coordinate_sd] = predict(centr_rss_gpr_y_model, per_tp_RSS);
        CENTR_RSS_GPR_PRED = [x_coordinate_tp y_coordinate_tp];
        dist_coords = [CENTR_RSS_GPR_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_rss(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        per_setup_err_ellipse_area_centr_rss(setup_idx,TP_idx) = 5.991*pi*x_coordinate_sd*y_coordinate_sd;
        
        a = x_coordinate_sd;
		b = y_coordinate_sd;
		x_diff = (TRUE_TP_COORDS(1) - x_coordinate_tp)^2 / a^2;
		y_diff = (TRUE_TP_COORDS(2) - y_coordinate_tp)^2 / b^2;
		
		if ((x_diff + y_diff) <= 5.991)
            NUM_TP_INSIDE_ELLIPSE_CR = NUM_TP_INSIDE_ELLIPSE_CR+1;
        end

        [x_coordinate_tp,x_coordinate_sd] = predict(centr_aoa_gpr_x_model, per_tp_est_angles);
        [y_coordinate_tp,y_coordinate_sd] = predict(centr_aoa_gpr_y_model, per_tp_est_angles);
        CENTR_AOA_GPR_PRED = [x_coordinate_tp y_coordinate_tp];
        dist_coords = [CENTR_AOA_GPR_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_aoa(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        per_setup_err_ellipse_area_centr_aoa(setup_idx,TP_idx) = 5.991*pi*x_coordinate_sd*y_coordinate_sd;
        
        a = x_coordinate_sd;
		b = y_coordinate_sd;
		x_diff = (TRUE_TP_COORDS(1) - x_coordinate_tp)^2 / a^2;
		y_diff = (TRUE_TP_COORDS(2) - y_coordinate_tp)^2 / b^2;
		
		if ((x_diff + y_diff) <= 5.991)
            NUM_TP_INSIDE_ELLIPSE_CA = NUM_TP_INSIDE_ELLIPSE_CA+1;
        end

        [x_coordinate_tp,x_coordinate_sd] = predict(centr_hybrid_gpr_x_model, [per_tp_RSS per_tp_est_angles]);
        [y_coordinate_tp,y_coordinate_sd] = predict(centr_hybrid_gpr_y_model, [per_tp_RSS per_tp_est_angles]);
        CENTR_HYBRID_GPR_PRED = [x_coordinate_tp y_coordinate_tp];
        dist_coords = [CENTR_HYBRID_GPR_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_hybrid(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        per_setup_err_ellipse_area_centr_hybrid(setup_idx,TP_idx) = 5.991*pi*x_coordinate_sd*y_coordinate_sd;

        a = x_coordinate_sd;
		b = y_coordinate_sd;
		x_diff = (TRUE_TP_COORDS(1) - x_coordinate_tp)^2 / a^2;
		y_diff = (TRUE_TP_COORDS(2) - y_coordinate_tp)^2 / b^2;
		
		if ((x_diff + y_diff) <= 5.991)
            NUM_TP_INSIDE_ELLIPSE_CH = NUM_TP_INSIDE_ELLIPSE_CH+1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%% DISTRIBUTED %%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %--------------- DISTRIBUTED MEDIAN ------------------------
        x_coordinate_median_tp = median(DISTR_PRED_TP_X_COORDS);
        y_coordinate_median_tp = median(DISTR_PRED_TP_Y_COORDS);
        DISTR_GPR_MEDIAN_PRED = [x_coordinate_median_tp y_coordinate_median_tp];
        dist_coords = [DISTR_GPR_MEDIAN_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_median(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');

        [sorted_X_COORDS, sorted_x_indices] = sort(DISTR_PRED_TP_X_COORDS);
        [sorted_Y_COORDS, sorted_y_indices] = sort(DISTR_PRED_TP_Y_COORDS);
        sorted_X_COORDS_SD = DISTR_PRED_TP_X_COORDS_SD(sorted_x_indices);
        sorted_Y_COORDS_SD = DISTR_PRED_TP_Y_COORDS_SD(sorted_y_indices);

        if mod(L, 2) == 1
            % If L is odd, find the median index
            median_index = (L + 1) / 2;
            x_coordinate_median_sd = sorted_X_COORDS_SD(median_index);
            y_coordinate_median_sd = sorted_Y_COORDS_SD(median_index);
        else
            % If L is even, find the two median indices
            median_index1 = L / 2;
            median_index2 = median_index1 + 1;
            median_sds_x = [sorted_X_COORDS_SD(median_index1), sorted_X_COORDS_SD(median_index2)];
            median_sds_y = [sorted_Y_COORDS_SD(median_index1), sorted_Y_COORDS_SD(median_index2)];
            x_coordinate_median_sd = 0.5*sqrt(sum(median_sds_x.^2));
            y_coordinate_median_sd = 0.5*sqrt(sum(median_sds_y.^2));
        end
        per_setup_err_ellipse_area_distr_median(setup_idx,TP_idx) = 5.991*pi*x_coordinate_median_sd*y_coordinate_median_sd;

        a = x_coordinate_median_sd;
		b = y_coordinate_median_sd;
		x_diff = (TRUE_TP_COORDS(1) - x_coordinate_median_tp)^2 / a^2;
		y_diff = (TRUE_TP_COORDS(2) - y_coordinate_median_tp)^2 / b^2;
		
		if ((x_diff + y_diff) <= 5.991)
            NUM_TP_INSIDE_ELLIPSE_DD = NUM_TP_INSIDE_ELLIPSE_DD + 1;
        end

        %--------------- DISTRIBUTED MEAN ------------------------
        x_coordinate_mean_tp = mean(DISTR_PRED_TP_X_COORDS);
        y_coordinate_mean_tp = mean(DISTR_PRED_TP_Y_COORDS);
        DISTR_GPR_MEAN_PRED = [x_coordinate_mean_tp y_coordinate_mean_tp];
        dist_coords = [DISTR_GPR_MEAN_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_mean(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');

        x_coordinate_mean_sd = sqrt(sum(DISTR_PRED_TP_X_COORDS_SD.^2))/L;
        y_coordinate_mean_sd = sqrt(sum(DISTR_PRED_TP_X_COORDS_SD.^2))/L;
        per_setup_err_ellipse_area_distr_mean(setup_idx,TP_idx) = 5.991*pi*x_coordinate_mean_sd*y_coordinate_mean_sd;

        a = x_coordinate_mean_sd;
		b = y_coordinate_mean_sd;
		x_diff = (TRUE_TP_COORDS(1) - x_coordinate_mean_tp)^2 / a^2;
		y_diff = (TRUE_TP_COORDS(2) - y_coordinate_mean_tp)^2 / b^2;
		
		if ((x_diff + y_diff) <= 5.991)
            NUM_TP_INSIDE_ELLIPSE_DM = NUM_TP_INSIDE_ELLIPSE_DM + 1;
        end
        %--------------- DISTRIBUTED MEAN Z-SCORE ------------------------
        x_coordinate_zscores = zscore(DISTR_PRED_TP_X_COORDS);
        y_coordinate_zscores = zscore(DISTR_PRED_TP_Y_COORDS);

        x_coordinate_zscore_mean_tp = mean(DISTR_PRED_TP_X_COORDS(abs(x_coordinate_zscores)<=0.6));
        y_coordinate_zscore_mean_tp = mean(DISTR_PRED_TP_Y_COORDS(abs(y_coordinate_zscores)<=0.6));
        DISTR_ZSCORE_GPR_PRED = [x_coordinate_zscore_mean_tp y_coordinate_zscore_mean_tp];
        dist_coords = [DISTR_ZSCORE_GPR_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_zscore_tz_0p6(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        
        % Find the indices where the absolute z-score is <= 0.6
        x_valid_indices = abs(x_coordinate_zscores) <= 0.6;
        y_valid_indices = abs(y_coordinate_zscores) <= 0.6;
        
        % Extract the corresponding SD values using logical indexing
        corresponding_sds_x = DISTR_PRED_TP_X_COORDS_SD(x_valid_indices);
        corresponding_sds_y = DISTR_PRED_TP_Y_COORDS_SD(y_valid_indices);
        x_coordinate_zscore_mean_sd = sqrt(sum(corresponding_sds_x.^2))/length(corresponding_sds_x);
        y_coordinate_zscore_mean_sd = sqrt(sum(corresponding_sds_y.^2))/length(corresponding_sds_y);
        per_setup_err_ellipse_area_distr_zscore_tz_0p6(setup_idx,TP_idx) = 5.991*pi*x_coordinate_zscore_mean_sd*y_coordinate_zscore_mean_sd;

        %----------------------------------------------------------------------------
        x_coordinate_zscore_mean_tp = mean(DISTR_PRED_TP_X_COORDS(abs(x_coordinate_zscores)<=0.8));
        y_coordinate_zscore_mean_tp = mean(DISTR_PRED_TP_Y_COORDS(abs(y_coordinate_zscores)<=0.8));
        DISTR_ZSCORE_GPR_PRED = [x_coordinate_zscore_mean_tp y_coordinate_zscore_mean_tp];
        dist_coords = [DISTR_ZSCORE_GPR_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_zscore_tz_0p8(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        
        % Find the indices where the absolute z-score is <= 0.8
        x_valid_indices = abs(x_coordinate_zscores) <= 0.8;
        y_valid_indices = abs(y_coordinate_zscores) <= 0.8;
        
        % Extract the corresponding SD values using logical indexing
        corresponding_sds_x = DISTR_PRED_TP_X_COORDS_SD(x_valid_indices);
        corresponding_sds_y = DISTR_PRED_TP_Y_COORDS_SD(y_valid_indices);
        x_coordinate_zscore_mean_sd = sqrt(sum(corresponding_sds_x.^2))/length(corresponding_sds_x);
        y_coordinate_zscore_mean_sd = sqrt(sum(corresponding_sds_y.^2))/length(corresponding_sds_y);
        per_setup_err_ellipse_area_distr_zscore_tz_0p8(setup_idx,TP_idx) = 5.991*pi*x_coordinate_zscore_mean_sd*y_coordinate_zscore_mean_sd;

        %----------------------------------------------------------------------------

        x_coordinate_zscore_mean_tp = mean(DISTR_PRED_TP_X_COORDS(abs(x_coordinate_zscores)<=1));
        y_coordinate_zscore_mean_tp = mean(DISTR_PRED_TP_Y_COORDS(abs(y_coordinate_zscores)<=1));
        DISTR_ZSCORE_GPR_PRED = [x_coordinate_zscore_mean_tp y_coordinate_zscore_mean_tp];
        dist_coords = [DISTR_ZSCORE_GPR_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_zscore_tz_1(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');

        % Find the indices where the absolute z-score is <= 1
        x_valid_indices = abs(x_coordinate_zscores) <= 1;
        y_valid_indices = abs(y_coordinate_zscores) <= 1;
        
        % Extract the corresponding SD values using logical indexing
        corresponding_sds_x = DISTR_PRED_TP_X_COORDS_SD(x_valid_indices);
        corresponding_sds_y = DISTR_PRED_TP_Y_COORDS_SD(y_valid_indices);
        x_coordinate_zscore_mean_sd = sqrt(sum(corresponding_sds_x.^2))/length(corresponding_sds_x);
        y_coordinate_zscore_mean_sd = sqrt(sum(corresponding_sds_y.^2))/length(corresponding_sds_y);
        per_setup_err_ellipse_area_distr_zscore_tz_1(setup_idx,TP_idx) = 5.991*pi*x_coordinate_zscore_mean_sd*y_coordinate_zscore_mean_sd;

        a = x_coordinate_zscore_mean_sd;
		b = y_coordinate_zscore_mean_sd;
		x_diff = (TRUE_TP_COORDS(1) - x_coordinate_zscore_mean_tp)^2 / a^2;
		y_diff = (TRUE_TP_COORDS(2) - y_coordinate_zscore_mean_tp)^2 / b^2;
		
		if ((x_diff + y_diff) <= 5.991)
            NUM_TP_INSIDE_ELLIPSE_DZ = NUM_TP_INSIDE_ELLIPSE_DZ + 1;
        end

        %----------------------------------------------------------------------------

        x_coordinate_zscore_mean_tp = mean(DISTR_PRED_TP_X_COORDS(abs(x_coordinate_zscores)<=1.2));
        y_coordinate_zscore_mean_tp = mean(DISTR_PRED_TP_Y_COORDS(abs(y_coordinate_zscores)<=1.2));
        DISTR_ZSCORE_GPR_PRED = [x_coordinate_zscore_mean_tp y_coordinate_zscore_mean_tp];
        dist_coords = [DISTR_ZSCORE_GPR_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_zscore_tz_1p2(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');

        % Find the indices where the absolute z-score is <= 1.2
        x_valid_indices = abs(x_coordinate_zscores) <= 1.2;
        y_valid_indices = abs(y_coordinate_zscores) <= 1.2;
        
        % Extract the corresponding SD values using logical indexing
        corresponding_sds_x = DISTR_PRED_TP_X_COORDS_SD(x_valid_indices);
        corresponding_sds_y = DISTR_PRED_TP_Y_COORDS_SD(y_valid_indices);
        x_coordinate_zscore_mean_sd = sqrt(sum(corresponding_sds_x.^2))/length(corresponding_sds_x);
        y_coordinate_zscore_mean_sd = sqrt(sum(corresponding_sds_y.^2))/length(corresponding_sds_y);
        per_setup_err_ellipse_area_distr_zscore_tz_1p2(setup_idx,TP_idx) = 5.991*pi*x_coordinate_zscore_mean_sd*y_coordinate_zscore_mean_sd;

        %----------------------------------------------------------------------------

        x_coordinate_zscore_mean_tp = mean(DISTR_PRED_TP_X_COORDS(abs(x_coordinate_zscores)<=1.4));
        y_coordinate_zscore_mean_tp = mean(DISTR_PRED_TP_Y_COORDS(abs(y_coordinate_zscores)<=1.4));
        DISTR_ZSCORE_GPR_PRED = [x_coordinate_zscore_mean_tp y_coordinate_zscore_mean_tp];
        dist_coords = [DISTR_ZSCORE_GPR_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_zscore_tz_1p4(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        
        % Find the indices where the absolute z-score is <= 1.4
        x_valid_indices = abs(x_coordinate_zscores) <= 1.4;
        y_valid_indices = abs(y_coordinate_zscores) <= 1.4;
        
        % Extract the corresponding SD values using logical indexing
        corresponding_sds_x = DISTR_PRED_TP_X_COORDS_SD(x_valid_indices);
        corresponding_sds_y = DISTR_PRED_TP_Y_COORDS_SD(y_valid_indices);
        x_coordinate_zscore_mean_sd = sqrt(sum(corresponding_sds_x.^2))/length(corresponding_sds_x);
        y_coordinate_zscore_mean_sd = sqrt(sum(corresponding_sds_y.^2))/length(corresponding_sds_y);
        per_setup_err_ellipse_area_distr_zscore_tz_1p4(setup_idx,TP_idx) = 5.991*pi*x_coordinate_zscore_mean_sd*y_coordinate_zscore_mean_sd;
        
        %--------------- DISTRIBUTED BAYESIAN ------------------------
        x_coordinate_bayesian_var = 1 / sum(1 ./ (DISTR_PRED_TP_X_COORDS_SD .^ 2));
        x_coordinate_bayesian_sd = sqrt(x_coordinate_bayesian_var);
        x_coordinate_bayesian_tp = x_coordinate_bayesian_var * sum(DISTR_PRED_TP_X_COORDS ./ (DISTR_PRED_TP_X_COORDS_SD .^ 2));

        y_coordinate_bayesian_var = 1 / sum(1 ./ (DISTR_PRED_TP_Y_COORDS_SD .^ 2));
        y_coordinate_bayesian_sd = sqrt(y_coordinate_bayesian_var);
        y_coordinate_bayesian_tp = y_coordinate_bayesian_var * sum(DISTR_PRED_TP_Y_COORDS ./ (DISTR_PRED_TP_Y_COORDS_SD .^ 2));

        DISTR_BAY_GPR_PRED = [x_coordinate_bayesian_tp y_coordinate_bayesian_tp];
        dist_coords = [DISTR_BAY_GPR_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_bayesian(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        per_setup_err_ellipse_area_distr_bayesian(setup_idx,TP_idx) = 5.991*pi*x_coordinate_bayesian_sd*y_coordinate_bayesian_sd;
        
        a = x_coordinate_bayesian_sd;
		b = y_coordinate_bayesian_sd;
		x_diff = (TRUE_TP_COORDS(1) - x_coordinate_bayesian_tp)^2 / a^2;
		y_diff = (TRUE_TP_COORDS(2) - y_coordinate_bayesian_tp)^2 / b^2;
		
		if ((x_diff + y_diff) <= 5.991)
            NUM_TP_INSIDE_ELLIPSE_DBAY = NUM_TP_INSIDE_ELLIPSE_DBAY + 1;
        end
        
        %--------------- DISTRIBUTED BAYESIAN Z-SCORE ------------------------
        % Find the indices where the absolute z-score is <= 0.6
        x_valid_indices = abs(x_coordinate_zscores) <= 0.6;
        y_valid_indices = abs(y_coordinate_zscores) <= 0.6;

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
        per_setup_positioning_err_distr_zscore_tz_0p6_bayesian(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        per_setup_err_ellipse_area_distr_zscore_tz_0p6_bayesian(setup_idx,TP_idx) = 5.991*pi*x_coordinate_zscore_bayesian_sd*y_coordinate_zscore_bayesian_sd;
        
        %----------------------------------------------------------------------------
        % Find the indices where the absolute z-score is <= 0.8
        x_valid_indices = abs(x_coordinate_zscores) <= 0.8;
        y_valid_indices = abs(y_coordinate_zscores) <= 0.8;

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
        per_setup_positioning_err_distr_zscore_tz_0p8_bayesian(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        per_setup_err_ellipse_area_distr_zscore_tz_0p8_bayesian(setup_idx,TP_idx) = 5.991*pi*x_coordinate_zscore_bayesian_sd*y_coordinate_zscore_bayesian_sd;
        %----------------------------------------------------------------------------

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
        per_setup_positioning_err_distr_zscore_tz_1_bayesian(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        per_setup_err_ellipse_area_distr_zscore_tz_1_bayesian(setup_idx,TP_idx) = 5.991*pi*x_coordinate_zscore_bayesian_sd*y_coordinate_zscore_bayesian_sd;
        
        a = x_coordinate_zscore_bayesian_sd;
		b = y_coordinate_zscore_bayesian_sd;
		x_diff = (TRUE_TP_COORDS(1) - x_coordinate_zscore_bayesian_tp)^2 / a^2;
		y_diff = (TRUE_TP_COORDS(2) - y_coordinate_zscore_bayesian_tp)^2 / b^2;
		
		if ((x_diff + y_diff) <= 5.991)
            NUM_TP_INSIDE_ELLIPSE_DZ_BAY = NUM_TP_INSIDE_ELLIPSE_DZ_BAY + 1;
        end
        %----------------------------------------------------------------------------

        % Find the indices where the absolute z-score is <= 1.2
        x_valid_indices = abs(x_coordinate_zscores) <= 1.2;
        y_valid_indices = abs(y_coordinate_zscores) <= 1.2;

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
        per_setup_positioning_err_distr_zscore_tz_1p2_bayesian(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        per_setup_err_ellipse_area_distr_zscore_tz_1p2_bayesian(setup_idx,TP_idx) = 5.991*pi*x_coordinate_zscore_bayesian_sd*y_coordinate_zscore_bayesian_sd;
        
        %----------------------------------------------------------------------------
        % Find the indices where the absolute z-score is <= 1.4
        x_valid_indices = abs(x_coordinate_zscores) <= 1.4;
        y_valid_indices = abs(y_coordinate_zscores) <= 1.4;

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
        per_setup_positioning_err_distr_zscore_tz_1p4_bayesian(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        per_setup_err_ellipse_area_distr_zscore_tz_1p4_bayesian(setup_idx,TP_idx) = 5.991*pi*x_coordinate_zscore_bayesian_sd*y_coordinate_zscore_bayesian_sd;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% CRLB %%%%%%%%%%%%%%%%%%%%%%%%%%%

        %--------------- CENTRALIZED AOA ------------------------
        [x_coordinate_tp_crlb,x_coordinate_sd_crlb] = predict(centr_aoa_gpr_x_model, per_tp_est_angles_crlb);
        [y_coordinate_tp_crlb,y_coordinate_sd_crlb] = predict(centr_aoa_gpr_y_model, per_tp_est_angles_crlb);
        CENTR_AOA_GPR_PRED_CRLB = [x_coordinate_tp_crlb y_coordinate_tp_crlb];
        dist_coords = [CENTR_AOA_GPR_PRED_CRLB;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_aoa_crlb(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        per_setup_err_ellipse_area_centr_aoa_crlb(setup_idx,TP_idx) = 5.991*pi*x_coordinate_sd_crlb*y_coordinate_sd_crlb;

        a = x_coordinate_sd_crlb;
		b = y_coordinate_sd_crlb;
		x_diff = (TRUE_TP_COORDS(1) - x_coordinate_tp_crlb)^2 / a^2;
		y_diff = (TRUE_TP_COORDS(2) - y_coordinate_tp_crlb)^2 / b^2;
		
		if ((x_diff + y_diff) <= 5.991)
            NUM_TP_INSIDE_ELLIPSE_CA_CRLB = NUM_TP_INSIDE_ELLIPSE_CA_CRLB + 1;
        end
       
        %--------------- CENTRALIZED HYBRID ------------------------

        [x_coordinate_tp_crlb,x_coordinate_sd_crlb] = predict(centr_hybrid_gpr_x_model, [per_tp_RSS per_tp_est_angles_crlb]);
        [y_coordinate_tp_crlb,y_coordinate_sd_crlb] = predict(centr_hybrid_gpr_y_model, [per_tp_RSS per_tp_est_angles_crlb]);
        CENTR_HYBRID_GPR_PRED_CRLB = [x_coordinate_tp_crlb y_coordinate_tp_crlb];
        dist_coords = [CENTR_HYBRID_GPR_PRED_CRLB;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_hybrid_crlb(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        per_setup_err_ellipse_area_centr_hybrid_crlb(setup_idx,TP_idx) = 5.991*pi*x_coordinate_sd_crlb*x_coordinate_sd_crlb;

        a = x_coordinate_sd_crlb;
		b = x_coordinate_sd_crlb;
		x_diff = (TRUE_TP_COORDS(1) - x_coordinate_tp_crlb)^2 / a^2;
		y_diff = (TRUE_TP_COORDS(2) - y_coordinate_tp_crlb)^2 / b^2;
		
		if ((x_diff + y_diff) <= 5.991)
            NUM_TP_INSIDE_ELLIPSE_CH_CRLB = NUM_TP_INSIDE_ELLIPSE_CH_CRLB + 1;
        end

        %--------------- DISTRIBUTED -----------------------------
        DISTR_PRED_TP_X_COORDS_CRLB = zeros(1,L); 
        DISTR_PRED_TP_Y_COORDS_CRLB = zeros(1,L); 
        DISTR_PRED_TP_X_COORDS_SD_CRLB = zeros(1,L); 
        DISTR_PRED_TP_Y_COORDS_SD_CRLB = zeros(1,L);
        for AP_idx = 1:L
            dyn_var_name_x = strcat('per_ap_gprx',num2str(AP_idx));
            dyn_var_name_y = strcat('per_ap_gpry',num2str(AP_idx));
            [x_coordinate_tp_crlb,x_coordinate_sd_crlb] = predict(variable.(dyn_var_name_x), [per_tp_RSS(AP_idx) per_tp_est_angles_crlb(AP_idx)]);
            [y_coordinate_tp_crlb,y_coordinate_sd_crlb] = predict(variable.(dyn_var_name_y), [per_tp_RSS(AP_idx) per_tp_est_angles_crlb(AP_idx)]);
            DISTR_PRED_TP_X_COORDS_CRLB(AP_idx) = x_coordinate_tp_crlb;
            DISTR_PRED_TP_Y_COORDS_CRLB(AP_idx) = y_coordinate_tp_crlb; 
            DISTR_PRED_TP_X_COORDS_SD_CRLB(AP_idx) = x_coordinate_sd_crlb;
            DISTR_PRED_TP_Y_COORDS_SD_CRLB(AP_idx) = y_coordinate_sd_crlb; 
        end

        %--------------- DISTRIBUTED MEDIAN ------------------------
        x_coordinate_median_tp_crlb = median(DISTR_PRED_TP_X_COORDS_CRLB);
        y_coordinate_median_tp_crlb = median(DISTR_PRED_TP_Y_COORDS_CRLB);
        DISTR_GPR_MEDIAN_PRED_CRLB = [x_coordinate_median_tp_crlb y_coordinate_median_tp_crlb];
        dist_coords = [DISTR_GPR_MEDIAN_PRED_CRLB;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_median_crlb(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');

        [sorted_X_COORDS_CRLB, sorted_x_indices_crlb] = sort(DISTR_PRED_TP_X_COORDS_CRLB);
        [sorted_Y_COORDS_CRLB, sorted_y_indices_crlb] = sort(DISTR_PRED_TP_Y_COORDS_CRLB);
        sorted_X_COORDS_SD_CRLB = DISTR_PRED_TP_X_COORDS_SD_CRLB(sorted_x_indices_crlb);
        sorted_Y_COORDS_SD_CRLB = DISTR_PRED_TP_Y_COORDS_SD_CRLB(sorted_y_indices_crlb);

        if mod(L, 2) == 1
            % If L is odd, find the median index
            median_index_crlb = (L + 1) / 2;
            x_coordinate_median_sd_crlb = sorted_X_COORDS_SD_CRLB(median_index_crlb);
            y_coordinate_median_sd_crlb = sorted_Y_COORDS_SD_CRLB(median_index_crlb);
        else
            % If L is even, find the two median indices
            median_index1_crlb = L / 2;
            median_index2_crlb = median_index1_crlb + 1;
            median_sds_x_crlb = [sorted_X_COORDS_SD_CRLB(median_index1_crlb), sorted_X_COORDS_SD_CRLB(median_index2_crlb)];
            median_sds_y_crlb = [sorted_Y_COORDS_SD_CRLB(median_index1_crlb), sorted_Y_COORDS_SD_CRLB(median_index2_crlb)];
            x_coordinate_median_sd_crlb = 0.5*sqrt(sum(median_sds_x_crlb.^2));
            y_coordinate_median_sd_crlb = 0.5*sqrt(sum(median_sds_y_crlb.^2));
        end
        per_setup_err_ellipse_area_distr_median_crlb(setup_idx,TP_idx) = 5.991*pi*x_coordinate_median_sd_crlb*y_coordinate_median_sd_crlb;

        a = x_coordinate_median_sd_crlb;
		b = y_coordinate_median_sd_crlb;
		x_diff = (TRUE_TP_COORDS(1) - x_coordinate_median_tp_crlb)^2 / a^2;
		y_diff = (TRUE_TP_COORDS(2) - y_coordinate_median_tp_crlb)^2 / b^2;
		
		if ((x_diff + y_diff) <= 5.991)
            NUM_TP_INSIDE_ELLIPSE_DD_CRLB = NUM_TP_INSIDE_ELLIPSE_DD_CRLB + 1;
        end

        %--------------- DISTRIBUTED MEAN ------------------------
        x_coordinate_mean_tp_crlb = mean(DISTR_PRED_TP_X_COORDS_CRLB);
        y_coordinate_mean_tp_crlb = mean(DISTR_PRED_TP_Y_COORDS_CRLB);
        DISTR_GPR_MEAN_PRED_CRLB = [x_coordinate_mean_tp_crlb y_coordinate_mean_tp_crlb];
        dist_coords = [DISTR_GPR_MEAN_PRED_CRLB;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_mean_crlb(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');

        x_coordinate_mean_sd_crlb = sqrt(sum(DISTR_PRED_TP_X_COORDS_SD_CRLB.^2))/L;
        y_coordinate_mean_sd_crlb = sqrt(sum(DISTR_PRED_TP_Y_COORDS_SD_CRLB.^2))/L;
        per_setup_err_ellipse_area_distr_mean_crlb(setup_idx,TP_idx) = 5.991*pi*x_coordinate_mean_sd_crlb*y_coordinate_mean_sd_crlb;

        a = x_coordinate_mean_sd_crlb;
		b = y_coordinate_mean_sd_crlb;
		x_diff = (TRUE_TP_COORDS(1) - x_coordinate_mean_tp_crlb)^2 / a^2;
		y_diff = (TRUE_TP_COORDS(2) - y_coordinate_mean_tp_crlb)^2 / b^2;
		
		if ((x_diff + y_diff) <= 5.991)
            NUM_TP_INSIDE_ELLIPSE_DM_CRLB = NUM_TP_INSIDE_ELLIPSE_DM_CRLB + 1;
        end
        %--------------- DISTRIBUTED MEAN Z-SCORE ------------------------
        x_coordinate_zscores_crlb = zscore(DISTR_PRED_TP_X_COORDS_CRLB);
        y_coordinate_zscores_crlb = zscore(DISTR_PRED_TP_Y_COORDS_CRLB);

        x_coordinate_zscore_mean_tp_crlb = mean(DISTR_PRED_TP_X_COORDS_CRLB(abs(x_coordinate_zscores_crlb)<=1));
        y_coordinate_zscore_mean_tp_crlb = mean(DISTR_PRED_TP_Y_COORDS_CRLB(abs(y_coordinate_zscores_crlb)<=1));
        DISTR_ZSCORE_GPR_PRED_CRLB = [x_coordinate_zscore_mean_tp_crlb y_coordinate_zscore_mean_tp_crlb];
        dist_coords = [DISTR_ZSCORE_GPR_PRED_CRLB;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_zscore_tz_1_crlb(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        
        % Find the indices where the absolute z-score is <= 1
        x_valid_indices_crlb = find(abs(x_coordinate_zscores_crlb) <= 1);
        y_valid_indices_crlb = find(abs(y_coordinate_zscores_crlb) <= 1);
        % Extract the corresponding SD values
        corresponding_sds_x_crlb = DISTR_PRED_TP_X_COORDS_SD_CRLB(x_valid_indices_crlb);
        corresponding_sds_y_crlb = DISTR_PRED_TP_Y_COORDS_SD_CRLB(y_valid_indices_crlb);
        x_coordinate_zscore_mean_sd_crlb = sqrt(sum(corresponding_sds_x_crlb.^2))/length(corresponding_sds_x_crlb);
        y_coordinate_zscore_mean_sd_crlb = sqrt(sum(corresponding_sds_y_crlb.^2))/length(corresponding_sds_y_crlb);
        per_setup_err_ellipse_area_distr_zscore_tz_1_crlb(setup_idx,TP_idx) = 5.991*pi*x_coordinate_zscore_mean_sd_crlb*y_coordinate_zscore_mean_sd_crlb;

        a = x_coordinate_zscore_mean_sd_crlb;
		b = y_coordinate_zscore_mean_sd_crlb;
		x_diff = (TRUE_TP_COORDS(1) - x_coordinate_zscore_mean_tp_crlb)^2 / a^2;
		y_diff = (TRUE_TP_COORDS(2) - y_coordinate_zscore_mean_tp_crlb)^2 / b^2;
		
		if ((x_diff + y_diff) <= 5.991)
            NUM_TP_INSIDE_ELLIPSE_DZ_CRLB = NUM_TP_INSIDE_ELLIPSE_DZ_CRLB + 1;
        end

        %--------------- DISTRIBUTED BAYESIAN ------------------------
        x_coordinate_bayesian_var_crlb = 1 / sum(1 ./ (DISTR_PRED_TP_X_COORDS_SD_CRLB .^ 2));
        x_coordinate_bayesian_sd_crlb = sqrt(x_coordinate_bayesian_var_crlb);
        x_coordinate_bayesian_tp_crlb = x_coordinate_bayesian_var_crlb * sum(DISTR_PRED_TP_X_COORDS_CRLB ./ (DISTR_PRED_TP_X_COORDS_SD_CRLB .^ 2));

        y_coordinate_bayesian_var_crlb = 1 / sum(1 ./ (DISTR_PRED_TP_Y_COORDS_SD_CRLB .^ 2));
        y_coordinate_bayesian_sd_crlb = sqrt(y_coordinate_bayesian_var_crlb);
        y_coordinate_bayesian_tp_crlb = y_coordinate_bayesian_var_crlb * sum(DISTR_PRED_TP_Y_COORDS_CRLB ./ (DISTR_PRED_TP_Y_COORDS_SD_CRLB .^ 2));

        DISTR_BAY_GPR_PRED_CRLB = [x_coordinate_bayesian_tp_crlb y_coordinate_bayesian_tp_crlb];
        dist_coords = [DISTR_BAY_GPR_PRED_CRLB;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_bayesian_crlb(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        per_setup_err_ellipse_area_distr_bayesian_crlb(setup_idx,TP_idx) = 5.991*pi*x_coordinate_bayesian_sd_crlb*y_coordinate_bayesian_sd_crlb;
        
        a = x_coordinate_bayesian_sd_crlb;
		b = y_coordinate_bayesian_sd_crlb;
		x_diff = (TRUE_TP_COORDS(1) - x_coordinate_bayesian_tp_crlb)^2 / a^2;
		y_diff = (TRUE_TP_COORDS(2) - y_coordinate_bayesian_tp_crlb)^2 / b^2;
		
		if ((x_diff + y_diff) <= 5.991)
            NUM_TP_INSIDE_ELLIPSE_DBAY_CRLB = NUM_TP_INSIDE_ELLIPSE_DBAY_CRLB + 1;
        end

        %--------------- DISTRIBUTED BAYESIAN Z-SCORE ------------------------
        % Find the indices where the absolute z-score is <= 1
        x_coordinate_zscores_crlb = zscore(DISTR_PRED_TP_X_COORDS_CRLB);
        y_coordinate_zscores_crlb = zscore(DISTR_PRED_TP_Y_COORDS_CRLB);
        x_valid_indices_crlb = abs(x_coordinate_zscores_crlb) <= 1;
        y_valid_indices_crlb = abs(y_coordinate_zscores_crlb) <= 1;

        % Extract the corresponding mean values using logical indexing
        corresponding_means_x_crlb = DISTR_PRED_TP_X_COORDS_CRLB(x_valid_indices_crlb);
        corresponding_means_y_crlb = DISTR_PRED_TP_Y_COORDS_CRLB(y_valid_indices_crlb);

        % Extract the corresponding SD values using logical indexing
        corresponding_sds_x_crlb = DISTR_PRED_TP_X_COORDS_SD_CRLB(x_valid_indices_crlb);
        corresponding_sds_y_crlb = DISTR_PRED_TP_Y_COORDS_SD(y_valid_indices_crlb);

        x_coordinate_zscore_bayesian_var_crlb = 1 / sum(1 ./ (corresponding_sds_x_crlb .^ 2));
        x_coordinate_zscore_bayesian_sd_crlb = sqrt(x_coordinate_zscore_bayesian_var_crlb);
        x_coordinate_zscore_bayesian_tp_crlb = x_coordinate_zscore_bayesian_var_crlb * sum(corresponding_means_x_crlb ./ (corresponding_sds_x_crlb .^ 2));

        y_coordinate_zscore_bayesian_var_crlb = 1 / sum(1 ./ (corresponding_sds_y_crlb .^ 2));
        y_coordinate_zscore_bayesian_sd_crlb = sqrt(y_coordinate_zscore_bayesian_var_crlb);
        y_coordinate_zscore_bayesian_tp_crlb = y_coordinate_zscore_bayesian_var_crlb * sum(corresponding_means_y_crlb ./ (corresponding_sds_y_crlb .^ 2));

        DISTR_BAY_ZSCORE_GPR_PRED_CRLB = [x_coordinate_zscore_bayesian_tp_crlb y_coordinate_zscore_bayesian_tp_crlb];
        dist_coords = [DISTR_BAY_ZSCORE_GPR_PRED_CRLB;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_zscore_tz_1_bayesian_crlb(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        per_setup_err_ellipse_area_distr_zscore_tz_1_bayesian_crlb(setup_idx,TP_idx) = 5.991*pi*x_coordinate_zscore_bayesian_sd_crlb*y_coordinate_zscore_bayesian_sd_crlb;
        
        a = x_coordinate_zscore_bayesian_sd_crlb;
		b = y_coordinate_zscore_bayesian_sd_crlb;
		x_diff = (TRUE_TP_COORDS(1) - x_coordinate_zscore_bayesian_tp_crlb)^2 / a^2;
		y_diff = (TRUE_TP_COORDS(2) - y_coordinate_zscore_bayesian_tp_crlb)^2 / b^2;
		
		if ((x_diff + y_diff) <= 5.991)
            NUM_TP_INSIDE_ELLIPSE_DZ_BAY_CRLB = NUM_TP_INSIDE_ELLIPSE_DZ_BAY_CRLB + 1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%% LR LOCALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%

        %--------------- CENTRALIZED ------------------------
        X_LR_test = [1 psi_fngprnt_tp(TP_idx,:)]; % Feature matrix for TPs
        y_pred_x = X_LR_test * w_x;
        y_pred_y = X_LR_test * w_y;
        LR_PRED = [y_pred_x y_pred_y];
        dist_coords = [LR_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_rss_LR(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');  
    
        X_test_LR_AOA = [1 nom_azi_angle_TP_music_est(TP_idx,:)]; % Feature matrix for TPs
        y_pred_x = X_test_LR_AOA * w_x_aoa;
        y_pred_y = X_test_LR_AOA * w_y_aoa;
        LR_PRED = [y_pred_x y_pred_y];
        dist_coords = [LR_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_aoa_LR(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');       
        
        X_test_LR_hybrid = [1 psi_fngprnt_tp(TP_idx,:) nom_azi_angle_TP_music_est(TP_idx,:)]; % Feature matrix for TPs
        y_pred_x = X_test_LR_hybrid * w_x_hybrid;
        y_pred_y = X_test_LR_hybrid * w_y_hybrid;
        LR_PRED = [y_pred_x y_pred_y];
        dist_coords = [LR_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_hybrid_LR(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');

        X_LR_test_dB = [1 psi_fngprnt_tp_dB(TP_idx,:)]; % Feature matrix for TPs
        y_pred_x = X_LR_test_dB * w_x_dB;
        y_pred_y = X_LR_test_dB * w_y_dB;
        LR_PRED = [y_pred_x y_pred_y];
        dist_coords = [LR_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_rss_LR_dB(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');  
    
        X_test_LR_hybrid_dB = [1 psi_fngprnt_tp_dB(TP_idx,:) nom_azi_angle_TP_music_est(TP_idx,:)]; % Feature matrix for TPs
        y_pred_x = X_test_LR_hybrid_dB * w_x_hybrid_dB;
        y_pred_y = X_test_LR_hybrid_dB * w_y_hybrid_dB;
        LR_PRED = [y_pred_x y_pred_y];
        dist_coords = [LR_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_hybrid_LR_dB(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');

        %--------------- DISTRIBUTED ------------------------
        DISTR_PRED_TP_X_COORDS_LR = zeros(1,L); 
        DISTR_PRED_TP_Y_COORDS_LR = zeros(1,L);    
        for AP_idx = 1:L
            X_test_LR_hybrid = [1 psi_fngprnt_tp_dB(TP_idx,AP_idx) nom_azi_angle_TP_music_est(TP_idx,AP_idx)];
            dyn_var_name_x = strcat('w_per_ap_x',num2str(AP_idx));
            dyn_var_name_y = strcat('w_per_ap_y',num2str(AP_idx));
            y_pred_x = X_test_LR_hybrid * variable.(dyn_var_name_x);
            y_pred_y = X_test_LR_hybrid * variable.(dyn_var_name_y);
            DISTR_PRED_TP_X_COORDS_LR(AP_idx) = y_pred_x;
            DISTR_PRED_TP_Y_COORDS_LR(AP_idx) = y_pred_y; 
        end

        y_pred_x_distr_median_tp = median(DISTR_PRED_TP_X_COORDS_LR);
        y_pred_y_distr_median_tp = median(DISTR_PRED_TP_Y_COORDS_LR);
        DISTR_MEDIAN_LR_PRED = [y_pred_x_distr_median_tp y_pred_y_distr_median_tp];
        dist_coords = [DISTR_MEDIAN_LR_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_median_LR_dB(setup_idx,TP_idx) = pdist(dist_coords,'euclidean'); 
            
        X_test_LR_AOA_crlb = [1 nom_azi_angle_TP_crlb(TP_idx,:)]; % Feature matrix for TPs
        y_pred_x = X_test_LR_AOA_crlb * w_x_aoa;
        y_pred_y = X_test_LR_AOA_crlb * w_y_aoa;
        LR_PRED_CRLB = [y_pred_x y_pred_y];
        dist_coords = [LR_PRED_CRLB;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_aoa_LR_crlb(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');       
        
        X_test_LR_hybrid_crlb = [1 psi_fngprnt_tp(TP_idx,:) nom_azi_angle_TP_crlb(TP_idx,:)]; % Feature matrix for TPs
        y_pred_x = X_test_LR_hybrid_crlb * w_x_hybrid;
        y_pred_y = X_test_LR_hybrid_crlb * w_y_hybrid;
        LR_PRED_CRLB = [y_pred_x y_pred_y];
        dist_coords = [LR_PRED_CRLB;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_hybrid_LR_crlb(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');

        X_test_LR_hybrid_dB_crlb = [1 psi_fngprnt_tp_dB(TP_idx,:) nom_azi_angle_TP_crlb(TP_idx,:)]; % Feature matrix for TPs
        y_pred_x = X_test_LR_hybrid_dB_crlb * w_x_hybrid_dB;
        y_pred_y = X_test_LR_hybrid_dB_crlb * w_y_hybrid_dB;
        LR_PRED_CRLB = [y_pred_x y_pred_y];
        dist_coords = [LR_PRED_CRLB;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_hybrid_LR_dB_crlb(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');

        DISTR_PRED_TP_X_COORDS_LR_CRLB = zeros(1,L); 
        DISTR_PRED_TP_Y_COORDS_LR_CRLB = zeros(1,L);    
        for AP_idx = 1:L
            X_test_LR_hybrid_crlb = [1 psi_fngprnt_tp_dB(TP_idx,AP_idx) nom_azi_angle_TP_crlb(TP_idx,AP_idx)];
            dyn_var_name_x = strcat('w_per_ap_x',num2str(AP_idx));
            dyn_var_name_y = strcat('w_per_ap_y',num2str(AP_idx));
            y_pred_x = X_test_LR_hybrid_crlb * variable.(dyn_var_name_x);
            y_pred_y = X_test_LR_hybrid_crlb * variable.(dyn_var_name_y);
            DISTR_PRED_TP_X_COORDS_LR_CRLB(AP_idx) = y_pred_x;
            DISTR_PRED_TP_Y_COORDS_LR_CRLB(AP_idx) = y_pred_y; 
        end

        y_pred_x_distr_median_tp_crlb = median(DISTR_PRED_TP_X_COORDS_LR_CRLB);
        y_pred_y_distr_median_tp_crlb = median(DISTR_PRED_TP_Y_COORDS_LR_CRLB);
        DISTR_MEDIAN_LR_PRED = [y_pred_x_distr_median_tp_crlb y_pred_y_distr_median_tp_crlb];
        dist_coords = [DISTR_MEDIAN_LR_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_median_LR_dB_crlb(setup_idx,TP_idx) = pdist(dist_coords,'euclidean'); 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% WKNN LOCALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %--------------- CENTRALIZED ------------------------
        X_test_KNN = psi_fngprnt_tp(TP_idx,:); % Feature matrix for TPs
        % Find the K nearest neighbors
        [idx, D] = knnsearch(X_train_KNN, X_test_KNN, 'K', KNN_K);
        % Calculate weighted average of nearest neighbors
        weights = 1 ./ D;
        weights = weights ./ sum(weights, 2); % Normalize weights

        y_pred = sum(y_train_KNN(idx) .* weights, 2);
        y_pred_x = real(y_pred);
        y_pred_y = imag(y_pred);
        KNN_PRED = [y_pred_x y_pred_y];
        dist_coords = [KNN_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_rss_WKNN(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');  

        X_test_KNN_AOA = nom_azi_angle_TP_music_est(TP_idx,:); % Feature matrix for TPs
        % Find the K nearest neighbors
        [idx, D] = knnsearch(X_train_KNN_AOA, X_test_KNN_AOA, 'K', KNN_K);
        % Calculate weighted average of nearest neighbors
        weights = 1 ./ D;
        weights = weights ./ sum(weights, 2); % Normalize weights
        
        y_pred = sum(y_train_KNN(idx) .* weights, 2);
        y_pred_x = real(y_pred);
        y_pred_y = imag(y_pred);
        KNN_PRED = [y_pred_x y_pred_y];
        dist_coords = [KNN_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_aoa_WKNN(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');  
        
        X_test_KNN_hybrid = [psi_fngprnt_tp(TP_idx,:) nom_azi_angle_TP_music_est(TP_idx,:)]; % Feature matrix for TPs
        %X_test_KNN_hybrid = [psi_fngprnt_tp(TP_idx,:)]; % Feature matrix for TPs
        % Find the K nearest neighbors
        [idx, D] = knnsearch(X_train_KNN_hybrid, X_test_KNN_hybrid, 'K', KNN_K);
        % Calculate weighted average of nearest neighbors
        weights = 1 ./ D;
        weights = weights ./ sum(weights, 2); % Normalize weights

        y_pred = sum(y_train_KNN(idx) .* weights, 2);
        y_pred_x = real(y_pred);
        y_pred_y = imag(y_pred);
        KNN_PRED = [y_pred_x y_pred_y];
        dist_coords = [KNN_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_hybrid_WKNN(setup_idx,TP_idx) = pdist(dist_coords,'euclidean'); 

        X_test_KNN_dB = psi_fngprnt_tp_dB(TP_idx,:); % Feature matrix for TPs
        % Find the K nearest neighbors
        [idx, D] = knnsearch(X_train_KNN_dB, X_test_KNN_dB, 'K', KNN_K);
        % Calculate weighted average of nearest neighbors
        weights = 1 ./ D;
        weights = weights ./ sum(weights, 2); % Normalize weights

        y_pred = sum(y_train_KNN(idx) .* weights, 2);
        y_pred_x = real(y_pred);
        y_pred_y = imag(y_pred);
        KNN_PRED = [y_pred_x y_pred_y];
        dist_coords = [KNN_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_rss_WKNN_dB(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');  

        X_test_KNN_hybrid_dB = [psi_fngprnt_tp_dB(TP_idx,:) nom_azi_angle_TP_music_est(TP_idx,:)]; % Feature matrix for TPs
        % Find the K nearest neighbors
        [idx, D] = knnsearch(X_train_KNN_hybrid_dB, X_test_KNN_hybrid_dB, 'K', KNN_K);
        % Calculate weighted average of nearest neighbors
        weights = 1 ./ D;
        weights = weights ./ sum(weights, 2); % Normalize weights

        y_pred = sum(y_train_KNN(idx) .* weights, 2);
        y_pred_x = real(y_pred);
        y_pred_y = imag(y_pred);
        KNN_PRED = [y_pred_x y_pred_y];
        dist_coords = [KNN_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_hybrid_WKNN_dB(setup_idx,TP_idx) = pdist(dist_coords,'euclidean'); 

        %--------------- DISTRIBUTED ------------------------
        DISTR_PRED_TP_X_COORDS_KNN = zeros(1,L); 
        DISTR_PRED_TP_Y_COORDS_KNN = zeros(1,L);
        for AP_idx = 1:L
            X_test_KNN_hybrid_per_ap = [psi_fngprnt_tp_dB(TP_idx,AP_idx) nom_azi_angle_TP_music_est(TP_idx,AP_idx)];
            dyn_var_name_knn = strcat('X_train_KNN_hybrid_per_ap',num2str(AP_idx));
            [idx, D] = knnsearch(variable.(dyn_var_name_knn), X_test_KNN_hybrid_per_ap, 'K', KNN_K);

            % Calculate weighted average of nearest neighbors
            weights = 1 ./ D;
            weights = weights ./ sum(weights, 2); % Normalize weights

            y_pred = sum(y_train_KNN(idx) .* weights, 2);
            y_pred_x = real(y_pred);
            y_pred_y = imag(y_pred);
            DISTR_PRED_TP_X_COORDS_KNN(AP_idx) = y_pred_x;
            DISTR_PRED_TP_Y_COORDS_KNN(AP_idx) = y_pred_y; 
        end
        y_pred_x_distr_median_tp = median(DISTR_PRED_TP_X_COORDS_KNN);
        y_pred_y_distr_median_tp = median(DISTR_PRED_TP_Y_COORDS_KNN);
        DISTR_MEDIAN_KNN_PRED = [y_pred_x_distr_median_tp y_pred_y_distr_median_tp];
        dist_coords = [DISTR_MEDIAN_KNN_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_median_WKNN_dB(setup_idx,TP_idx) = pdist(dist_coords,'euclidean'); 
            
        X_test_KNN_AOA_crlb = [nom_azi_angle_TP_crlb(TP_idx,:)]; % Feature matrix for TPs
        % Find the K nearest neighbors
        [idx, D] = knnsearch(X_train_KNN_AOA_crlb, X_test_KNN_AOA_crlb, 'K', KNN_K);
        % Calculate weighted average of nearest neighbors
        weights = 1 ./ D;
        weights = weights ./ sum(weights, 2); % Normalize weights
            
        y_pred = sum(y_train_KNN(idx) .* weights, 2);
        y_pred_x = real(y_pred);
        y_pred_y = imag(y_pred);
        KNN_PRED = [y_pred_x y_pred_y];
        dist_coords = [KNN_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_aoa_WKNN_crlb(setup_idx,TP_idx) = pdist(dist_coords,'euclidean'); 

        X_test_KNN_hybrid_dB_crlb = [psi_fngprnt_tp_dB(TP_idx,:) nom_azi_angle_TP_crlb(TP_idx,:)]; % Feature matrix for TPs
        % Find the K nearest neighbors
        [idx, D] = knnsearch(X_train_KNN_hybrid_dB_crlb, X_test_KNN_hybrid_dB_crlb, 'K', KNN_K);
        % Calculate weighted average of nearest neighbors
        weights = 1 ./ D;
        weights = weights ./ sum(weights, 2); % Normalize weights
            
        y_pred = sum(y_train_KNN(idx) .* weights, 2);
        y_pred_x = real(y_pred);
        y_pred_y = imag(y_pred);
        KNN_PRED = [y_pred_x y_pred_y];
        dist_coords = [KNN_PRED;TRUE_TP_COORDS];
        per_setup_positioning_err_centr_hybrid_WKNN_dB_crlb(setup_idx,TP_idx) = pdist(dist_coords,'euclidean');
        
        DISTR_PRED_TP_X_COORDS_KNN_CRLB = zeros(1,L); 
        DISTR_PRED_TP_Y_COORDS_KNN_CRLB = zeros(1,L);
        for AP_idx = 1:L
            X_test_KNN_hybrid_per_ap_crlb = [psi_fngprnt_tp_dB(TP_idx,AP_idx) nom_azi_angle_TP_crlb(TP_idx,AP_idx)];
            dyn_var_name_knn = strcat('X_train_KNN_hybrid_per_ap',num2str(AP_idx));
            [idx, D] = knnsearch(variable.(dyn_var_name_knn), X_test_KNN_hybrid_per_ap_crlb, 'K', KNN_K);

            % Calculate weighted average of nearest neighbors
            weights = 1 ./ D;
            weights = weights ./ sum(weights, 2); % Normalize weights

            y_pred = sum(y_train_KNN(idx) .* weights, 2);
            y_pred_x = real(y_pred);
            y_pred_y = imag(y_pred);
            DISTR_PRED_TP_X_COORDS_KNN_CRLB(AP_idx) = y_pred_x;
            DISTR_PRED_TP_Y_COORDS_KNN_CRLB(AP_idx) = y_pred_y; 
        end
        y_pred_x_distr_median_tp_crlb = median(DISTR_PRED_TP_X_COORDS_KNN_CRLB);
        y_pred_y_distr_median_tp_crlb = median(DISTR_PRED_TP_Y_COORDS_KNN_CRLB);
        DISTR_MEDIAN_KNN_PRED_CRLB = [y_pred_x_distr_median_tp_crlb y_pred_y_distr_median_tp_crlb];
        dist_coords = [DISTR_MEDIAN_KNN_PRED_CRLB;TRUE_TP_COORDS];
        per_setup_positioning_err_distr_median_WKNN_dB_crlb(setup_idx,TP_idx) = pdist(dist_coords,'euclidean'); 
            
    end %TP_idx = 1:num_tp_points 

    % Prepare test input for FCNN
    X_TP = nom_azi_angle_TP_music_est';  % Transpose to 25 x 1000
    
    % Predict (x, y) coordinates using trained network
    Y_TP_pred = fcnn_model(X_TP);               % Size: 2 x 1000
    x_TP_pred = Y_TP_pred(1, :)';
    y_TP_pred = Y_TP_pred(2, :)';
    
    % Ground truth coordinates from complex positions
    x_TP_true = real(TP_positions(:));   % Ensure column vector
    y_TP_true = imag(TP_positions(:));
    
    % Compute Euclidean error for each test point
    errors_TP = sqrt((x_TP_pred - x_TP_true).^2 + (y_TP_pred - y_TP_true).^2);
    
    % Compute mean error
    positioning_err_fcnn(setup_idx,1) = mean(errors_TP);
    
    end %setup_idx = 1:nbrOfSetups

% Average positioning errors for centralized approaches
avg_positioning_err_centr_rss(n_idx) = mean(per_setup_positioning_err_centr_rss(~isnan(per_setup_positioning_err_centr_rss) & ~isinf(per_setup_positioning_err_centr_rss)));
avg_positioning_err_centr_aoa(n_idx) = mean(per_setup_positioning_err_centr_aoa(~isnan(per_setup_positioning_err_centr_aoa) & ~isinf(per_setup_positioning_err_centr_aoa)));
avg_positioning_err_centr_hybrid(n_idx) = mean(per_setup_positioning_err_centr_hybrid(~isnan(per_setup_positioning_err_centr_hybrid) & ~isinf(per_setup_positioning_err_centr_hybrid)));

% Average positioning errors for distributed approaches
avg_positioning_err_distr_median(n_idx) = mean(per_setup_positioning_err_distr_median(~isnan(per_setup_positioning_err_distr_median) & ~isinf(per_setup_positioning_err_distr_median)));
avg_positioning_err_distr_mean(n_idx) = mean(per_setup_positioning_err_distr_mean(~isnan(per_setup_positioning_err_distr_mean) & ~isinf(per_setup_positioning_err_distr_mean)));
avg_positioning_err_distr_zscore_tz_0p6(n_idx) = mean(per_setup_positioning_err_distr_zscore_tz_0p6(~isnan(per_setup_positioning_err_distr_zscore_tz_0p6) & ~isinf(per_setup_positioning_err_distr_zscore_tz_0p6)));
avg_positioning_err_distr_zscore_tz_0p8(n_idx) = mean(per_setup_positioning_err_distr_zscore_tz_0p8(~isnan(per_setup_positioning_err_distr_zscore_tz_0p8) & ~isinf(per_setup_positioning_err_distr_zscore_tz_0p8)));
avg_positioning_err_distr_zscore_tz_1(n_idx) = mean(per_setup_positioning_err_distr_zscore_tz_1(~isnan(per_setup_positioning_err_distr_zscore_tz_1) & ~isinf(per_setup_positioning_err_distr_zscore_tz_1)));
avg_positioning_err_distr_zscore_tz_1p2(n_idx) = mean(per_setup_positioning_err_distr_zscore_tz_1p2(~isnan(per_setup_positioning_err_distr_zscore_tz_1p2) & ~isinf(per_setup_positioning_err_distr_zscore_tz_1p2)));
avg_positioning_err_distr_zscore_tz_1p4(n_idx) = mean(per_setup_positioning_err_distr_zscore_tz_1p4(~isnan(per_setup_positioning_err_distr_zscore_tz_1p4) & ~isinf(per_setup_positioning_err_distr_zscore_tz_1p4)));
avg_positioning_err_distr_bayesian(n_idx) = mean(per_setup_positioning_err_distr_bayesian(~isnan(per_setup_positioning_err_distr_bayesian) & ~isinf(per_setup_positioning_err_distr_bayesian)));
avg_positioning_err_distr_zscore_tz_0p6_bayesian(n_idx) = mean(per_setup_positioning_err_distr_zscore_tz_0p6_bayesian(~isnan(per_setup_positioning_err_distr_zscore_tz_0p6_bayesian) & ~isinf(per_setup_positioning_err_distr_zscore_tz_0p6_bayesian)));
avg_positioning_err_distr_zscore_tz_0p8_bayesian(n_idx) = mean(per_setup_positioning_err_distr_zscore_tz_0p8_bayesian(~isnan(per_setup_positioning_err_distr_zscore_tz_0p8_bayesian) & ~isinf(per_setup_positioning_err_distr_zscore_tz_0p8_bayesian)));
avg_positioning_err_distr_zscore_tz_1_bayesian(n_idx) = mean(per_setup_positioning_err_distr_zscore_tz_1_bayesian(~isnan(per_setup_positioning_err_distr_zscore_tz_1_bayesian) & ~isinf(per_setup_positioning_err_distr_zscore_tz_1_bayesian)));
avg_positioning_err_distr_zscore_tz_1p2_bayesian(n_idx) = mean(per_setup_positioning_err_distr_zscore_tz_1p2_bayesian(~isnan(per_setup_positioning_err_distr_zscore_tz_1p2_bayesian) & ~isinf(per_setup_positioning_err_distr_zscore_tz_1p2_bayesian)));
avg_positioning_err_distr_zscore_tz_1p4_bayesian(n_idx) = mean(per_setup_positioning_err_distr_zscore_tz_1p4_bayesian(~isnan(per_setup_positioning_err_distr_zscore_tz_1p4_bayesian) & ~isinf(per_setup_positioning_err_distr_zscore_tz_1p4_bayesian)));

% Average positioning errors for CRLB (Cramér-Rao lower bound) in various approaches
avg_positioning_err_centr_hybrid_crlb(n_idx) = mean(per_setup_positioning_err_centr_hybrid_crlb(~isnan(per_setup_positioning_err_centr_hybrid_crlb) & ~isinf(per_setup_positioning_err_centr_hybrid_crlb)));
avg_positioning_err_centr_aoa_crlb(n_idx) = mean(per_setup_positioning_err_centr_aoa_crlb(~isnan(per_setup_positioning_err_centr_aoa_crlb) & ~isinf(per_setup_positioning_err_centr_aoa_crlb)));
avg_positioning_err_distr_median_crlb(n_idx) = mean(per_setup_positioning_err_distr_median_crlb(~isnan(per_setup_positioning_err_distr_median_crlb) & ~isinf(per_setup_positioning_err_distr_median_crlb)));
avg_positioning_err_distr_mean_crlb(n_idx) = mean(per_setup_positioning_err_distr_mean_crlb(~isnan(per_setup_positioning_err_distr_mean_crlb) & ~isinf(per_setup_positioning_err_distr_mean_crlb)));
avg_positioning_err_distr_zscore_tz_1_crlb(n_idx) = mean(per_setup_positioning_err_distr_zscore_tz_1_crlb(~isnan(per_setup_positioning_err_distr_zscore_tz_1_crlb) & ~isinf(per_setup_positioning_err_distr_zscore_tz_1_crlb)));
avg_positioning_err_distr_bayesian_crlb(n_idx) = mean(per_setup_positioning_err_distr_bayesian_crlb(~isnan(per_setup_positioning_err_distr_bayesian_crlb) & ~isinf(per_setup_positioning_err_distr_bayesian_crlb)));
avg_positioning_err_distr_zscore_tz_1_bayesian_crlb(n_idx) = mean(per_setup_positioning_err_distr_zscore_tz_1_bayesian_crlb(~isnan(per_setup_positioning_err_distr_zscore_tz_1_bayesian_crlb) & ~isinf(per_setup_positioning_err_distr_zscore_tz_1_bayesian_crlb)));

% Average error ellipse area for centralized approaches
avg_err_ellipse_area_centr_rss(n_idx) = mean(per_setup_err_ellipse_area_centr_rss(~isnan(per_setup_err_ellipse_area_centr_rss) & ~isinf(per_setup_err_ellipse_area_centr_rss)));
avg_err_ellipse_area_centr_aoa(n_idx) = mean(per_setup_err_ellipse_area_centr_aoa(~isnan(per_setup_err_ellipse_area_centr_aoa) & ~isinf(per_setup_err_ellipse_area_centr_aoa)));
avg_err_ellipse_area_centr_hybrid(n_idx) = mean(per_setup_err_ellipse_area_centr_hybrid(~isnan(per_setup_err_ellipse_area_centr_hybrid) & ~isinf(per_setup_err_ellipse_area_centr_hybrid)));

% Average error ellipse area for distributed approaches
avg_err_ellipse_area_distr_median(n_idx) = mean(per_setup_err_ellipse_area_distr_median(~isnan(per_setup_err_ellipse_area_distr_median) & ~isinf(per_setup_err_ellipse_area_distr_median)));
avg_err_ellipse_area_distr_mean(n_idx) = mean(per_setup_err_ellipse_area_distr_mean(~isnan(per_setup_err_ellipse_area_distr_mean) & ~isinf(per_setup_err_ellipse_area_distr_mean)));
avg_err_ellipse_area_distr_zscore_tz_0p6(n_idx) = mean(per_setup_err_ellipse_area_distr_zscore_tz_0p6(~isnan(per_setup_err_ellipse_area_distr_zscore_tz_0p6) & ~isinf(per_setup_err_ellipse_area_distr_zscore_tz_0p6)));
avg_err_ellipse_area_distr_zscore_tz_0p8(n_idx) = mean(per_setup_err_ellipse_area_distr_zscore_tz_0p8(~isnan(per_setup_err_ellipse_area_distr_zscore_tz_0p8) & ~isinf(per_setup_err_ellipse_area_distr_zscore_tz_0p8)));
avg_err_ellipse_area_distr_zscore_tz_1(n_idx) = mean(per_setup_err_ellipse_area_distr_zscore_tz_1(~isnan(per_setup_err_ellipse_area_distr_zscore_tz_1) & ~isinf(per_setup_err_ellipse_area_distr_zscore_tz_1)));
avg_err_ellipse_area_distr_zscore_tz_1p2(n_idx) = mean(per_setup_err_ellipse_area_distr_zscore_tz_1p2(~isnan(per_setup_err_ellipse_area_distr_zscore_tz_1p2) & ~isinf(per_setup_err_ellipse_area_distr_zscore_tz_1p2)));
avg_err_ellipse_area_distr_zscore_tz_1p4(n_idx) = mean(per_setup_err_ellipse_area_distr_zscore_tz_1p4(~isnan(per_setup_err_ellipse_area_distr_zscore_tz_1p4) & ~isinf(per_setup_err_ellipse_area_distr_zscore_tz_1p4)));
avg_err_ellipse_area_distr_bayesian(n_idx) = mean(per_setup_err_ellipse_area_distr_bayesian(~isnan(per_setup_err_ellipse_area_distr_bayesian) & ~isinf(per_setup_err_ellipse_area_distr_bayesian)));
avg_err_ellipse_area_distr_zscore_tz_0p6_bayesian(n_idx) = mean(per_setup_err_ellipse_area_distr_zscore_tz_0p6_bayesian(~isnan(per_setup_err_ellipse_area_distr_zscore_tz_0p6_bayesian) & ~isinf(per_setup_err_ellipse_area_distr_zscore_tz_0p6_bayesian)));
avg_err_ellipse_area_distr_zscore_tz_0p8_bayesian(n_idx) = mean(per_setup_err_ellipse_area_distr_zscore_tz_0p8_bayesian(~isnan(per_setup_err_ellipse_area_distr_zscore_tz_0p8_bayesian) & ~isinf(per_setup_err_ellipse_area_distr_zscore_tz_0p8_bayesian)));
avg_err_ellipse_area_distr_zscore_tz_1_bayesian(n_idx) = mean(per_setup_err_ellipse_area_distr_zscore_tz_1_bayesian(~isnan(per_setup_err_ellipse_area_distr_zscore_tz_1_bayesian) & ~isinf(per_setup_err_ellipse_area_distr_zscore_tz_1_bayesian)));
avg_err_ellipse_area_distr_zscore_tz_1p2_bayesian(n_idx) = mean(per_setup_err_ellipse_area_distr_zscore_tz_1p2_bayesian(~isnan(per_setup_err_ellipse_area_distr_zscore_tz_1p2_bayesian) & ~isinf(per_setup_err_ellipse_area_distr_zscore_tz_1p2_bayesian)));
avg_err_ellipse_area_distr_zscore_tz_1p4_bayesian(n_idx) = mean(per_setup_err_ellipse_area_distr_zscore_tz_1p4_bayesian(~isnan(per_setup_err_ellipse_area_distr_zscore_tz_1p4_bayesian) & ~isinf(per_setup_err_ellipse_area_distr_zscore_tz_1p4_bayesian)));

% Average error ellipse area for CRLB in various approaches
avg_err_ellipse_area_centr_hybrid_crlb(n_idx) = mean(per_setup_err_ellipse_area_centr_hybrid_crlb(~isnan(per_setup_err_ellipse_area_centr_hybrid_crlb) & ~isinf(per_setup_err_ellipse_area_centr_hybrid_crlb)));
avg_err_ellipse_area_centr_aoa_crlb(n_idx) = mean(per_setup_err_ellipse_area_centr_aoa_crlb(~isnan(per_setup_err_ellipse_area_centr_aoa_crlb) & ~isinf(per_setup_err_ellipse_area_centr_aoa_crlb)));
avg_err_ellipse_area_distr_median_crlb(n_idx) = mean(per_setup_err_ellipse_area_distr_median_crlb(~isnan(per_setup_err_ellipse_area_distr_median_crlb) & ~isinf(per_setup_err_ellipse_area_distr_median_crlb)));
avg_err_ellipse_area_distr_mean_crlb(n_idx) = mean(per_setup_err_ellipse_area_distr_mean_crlb(~isnan(per_setup_err_ellipse_area_distr_mean_crlb) & ~isinf(per_setup_err_ellipse_area_distr_mean_crlb)));
avg_err_ellipse_area_distr_zscore_tz_1_crlb(n_idx) = mean(per_setup_err_ellipse_area_distr_zscore_tz_1_crlb(~isnan(per_setup_err_ellipse_area_distr_zscore_tz_1_crlb) & ~isinf(per_setup_err_ellipse_area_distr_zscore_tz_1_crlb)));
avg_err_ellipse_area_distr_bayesian_crlb(n_idx) = mean(per_setup_err_ellipse_area_distr_bayesian_crlb(~isnan(per_setup_err_ellipse_area_distr_bayesian_crlb) & ~isinf(per_setup_err_ellipse_area_distr_bayesian_crlb)));
avg_err_ellipse_area_distr_zscore_tz_1_bayesian_crlb(n_idx) = mean(per_setup_err_ellipse_area_distr_zscore_tz_1_bayesian_crlb(~isnan(per_setup_err_ellipse_area_distr_zscore_tz_1_bayesian_crlb) & ~isinf(per_setup_err_ellipse_area_distr_zscore_tz_1_bayesian_crlb)));

% Average positioning errors for LR (Linear Regression) approaches
avg_positioning_err_centr_rss_LR(n_idx) = mean(per_setup_positioning_err_centr_rss_LR(~isnan(per_setup_positioning_err_centr_rss_LR) & ~isinf(per_setup_positioning_err_centr_rss_LR)));
avg_positioning_err_centr_aoa_LR(n_idx) = mean(per_setup_positioning_err_centr_aoa_LR(~isnan(per_setup_positioning_err_centr_aoa_LR) & ~isinf(per_setup_positioning_err_centr_aoa_LR)));
avg_positioning_err_centr_hybrid_LR(n_idx) = mean(per_setup_positioning_err_centr_hybrid_LR(~isnan(per_setup_positioning_err_centr_hybrid_LR) & ~isinf(per_setup_positioning_err_centr_hybrid_LR)));

avg_positioning_err_centr_rss_LR_dB(n_idx) = mean(per_setup_positioning_err_centr_rss_LR_dB(~isnan(per_setup_positioning_err_centr_rss_LR_dB) & ~isinf(per_setup_positioning_err_centr_rss_LR_dB)));
avg_positioning_err_centr_hybrid_LR_dB(n_idx) = mean(per_setup_positioning_err_centr_hybrid_LR_dB(~isnan(per_setup_positioning_err_centr_hybrid_LR_dB) & ~isinf(per_setup_positioning_err_centr_hybrid_LR_dB)));
avg_positioning_err_distr_median_LR_dB(n_idx) = mean(per_setup_positioning_err_distr_median_LR_dB(~isnan(per_setup_positioning_err_distr_median_LR_dB) & ~isinf(per_setup_positioning_err_distr_median_LR_dB)));

% CRLB for LR approaches
avg_positioning_err_centr_aoa_LR_crlb(n_idx) = mean(per_setup_positioning_err_centr_aoa_LR_crlb(~isnan(per_setup_positioning_err_centr_aoa_LR_crlb) & ~isinf(per_setup_positioning_err_centr_aoa_LR_crlb)));
avg_positioning_err_centr_hybrid_LR_crlb(n_idx) = mean(per_setup_positioning_err_centr_hybrid_LR_crlb(~isnan(per_setup_positioning_err_centr_hybrid_LR_crlb) & ~isinf(per_setup_positioning_err_centr_hybrid_LR_crlb)));
avg_positioning_err_centr_hybrid_LR_dB_crlb(n_idx) = mean(per_setup_positioning_err_centr_hybrid_LR_dB_crlb(~isnan(per_setup_positioning_err_centr_hybrid_LR_dB_crlb) & ~isinf(per_setup_positioning_err_centr_hybrid_LR_dB_crlb)));
avg_positioning_err_distr_median_LR_dB_crlb(n_idx) = mean(per_setup_positioning_err_distr_median_LR_dB_crlb(~isnan(per_setup_positioning_err_distr_median_LR_dB_crlb) & ~isinf(per_setup_positioning_err_distr_median_LR_dB_crlb)));

% Corrected average positioning errors for WKNN approaches
avg_positioning_err_centr_rss_WKNN(n_idx) = mean(per_setup_positioning_err_centr_rss_WKNN(~isnan(per_setup_positioning_err_centr_rss_WKNN) & ~isinf(per_setup_positioning_err_centr_rss_WKNN)));
avg_positioning_err_centr_aoa_WKNN(n_idx) = mean(per_setup_positioning_err_centr_aoa_WKNN(~isnan(per_setup_positioning_err_centr_aoa_WKNN) & ~isinf(per_setup_positioning_err_centr_aoa_WKNN)));
avg_positioning_err_centr_hybrid_WKNN(n_idx) = mean(per_setup_positioning_err_centr_hybrid_WKNN(~isnan(per_setup_positioning_err_centr_hybrid_WKNN) & ~isinf(per_setup_positioning_err_centr_hybrid_WKNN)));

% Corrected average positioning errors for WKNN approaches in dB
avg_positioning_err_centr_rss_WKNN_dB(n_idx) = mean(per_setup_positioning_err_centr_rss_WKNN_dB(~isnan(per_setup_positioning_err_centr_rss_WKNN_dB) & ~isinf(per_setup_positioning_err_centr_rss_WKNN_dB)));
avg_positioning_err_centr_hybrid_WKNN_dB(n_idx) = mean(per_setup_positioning_err_centr_hybrid_WKNN_dB(~isnan(per_setup_positioning_err_centr_hybrid_WKNN_dB) & ~isinf(per_setup_positioning_err_centr_hybrid_WKNN_dB)));
avg_positioning_err_distr_median_WKNN_dB(n_idx) = mean(per_setup_positioning_err_distr_median_WKNN_dB(~isnan(per_setup_positioning_err_distr_median_WKNN_dB) & ~isinf(per_setup_positioning_err_distr_median_WKNN_dB)));

% Corrected CRLB for WKNN approaches
avg_positioning_err_centr_aoa_WKNN_crlb(n_idx) = mean(per_setup_positioning_err_centr_aoa_WKNN_crlb(~isnan(per_setup_positioning_err_centr_aoa_WKNN_crlb) & ~isinf(per_setup_positioning_err_centr_aoa_WKNN_crlb)));
avg_positioning_err_centr_hybrid_WKNN_dB_crlb(n_idx) = mean(per_setup_positioning_err_centr_hybrid_WKNN_dB_crlb(~isnan(per_setup_positioning_err_centr_hybrid_WKNN_dB_crlb) & ~isinf(per_setup_positioning_err_centr_hybrid_WKNN_dB_crlb)));
avg_positioning_err_distr_median_WKNN_dB_crlb(n_idx) = mean(per_setup_positioning_err_distr_median_WKNN_dB_crlb(~isnan(per_setup_positioning_err_distr_median_WKNN_dB_crlb) & ~isinf(per_setup_positioning_err_distr_median_WKNN_dB_crlb)));

% Average positioning error for FCNN
avg_positioning_err_fcnn(n_idx) = mean(positioning_err_fcnn);

% Print average positioning errors
fprintf('\nAverage Positioning Errors:\n');

% Centralized Positioning Errors
fprintf('C-RSS: %.4f\n', avg_positioning_err_centr_rss(n_idx));
fprintf('C-AOA: %.4f\n', avg_positioning_err_centr_aoa(n_idx));
fprintf('C-Hybrid: %.4f\n', avg_positioning_err_centr_hybrid(n_idx));

% Distributed Positioning Errors
fprintf('D-Median: %.4f\n', avg_positioning_err_distr_median(n_idx));
fprintf('D-Mean: %.4f\n', avg_positioning_err_distr_mean(n_idx));
fprintf('D-ZS TZ 0.6: %.4f\n', avg_positioning_err_distr_zscore_tz_0p6(n_idx));
fprintf('D-ZS TZ 0.8: %.4f\n', avg_positioning_err_distr_zscore_tz_0p8(n_idx));
fprintf('D-ZS TZ 1: %.4f\n', avg_positioning_err_distr_zscore_tz_1(n_idx));
fprintf('D-ZS TZ 1.2: %.4f\n', avg_positioning_err_distr_zscore_tz_1p2(n_idx));
fprintf('D-ZS TZ 1.4: %.4f\n', avg_positioning_err_distr_zscore_tz_1p4(n_idx));
fprintf('D-BAY: %.4f\n', avg_positioning_err_distr_bayesian(n_idx));
fprintf('D-ZS-BAY TZ 0.6: %.4f\n', avg_positioning_err_distr_zscore_tz_0p6_bayesian(n_idx));
fprintf('D-ZS-BAY TZ 0.8: %.4f\n', avg_positioning_err_distr_zscore_tz_0p8_bayesian(n_idx));
fprintf('D-ZS-BAY TZ 1: %.4f\n', avg_positioning_err_distr_zscore_tz_1_bayesian(n_idx));
fprintf('D-ZS-BAY TZ 1.2: %.4f\n', avg_positioning_err_distr_zscore_tz_1p2_bayesian(n_idx));
fprintf('D-ZS-BAY TZ 1.4: %.4f\n', avg_positioning_err_distr_zscore_tz_1p4_bayesian(n_idx));

% CRLB Errors
fprintf('C-Hybrid (CRLB): %.4f\n', avg_positioning_err_centr_hybrid_crlb(n_idx));
fprintf('C-AOA (CRLB): %.4f\n', avg_positioning_err_centr_aoa_crlb(n_idx));
fprintf('D-Median (CRLB): %.4f\n', avg_positioning_err_distr_median_crlb(n_idx));
fprintf('D-Mean (CRLB): %.4f\n', avg_positioning_err_distr_mean_crlb(n_idx));
fprintf('D-ZS TZ 1 (CRLB): %.4f\n', avg_positioning_err_distr_zscore_tz_1_crlb(n_idx));
fprintf('D-BAY (CRLB): %.4f\n', avg_positioning_err_distr_bayesian_crlb(n_idx));
fprintf('D-ZS-BAY TZ 1 (CRLB): %.4f\n', avg_positioning_err_distr_zscore_tz_1_bayesian_crlb(n_idx));

% Linear Regression (LR)
fprintf('C-RSS (LR): %.4f\n', avg_positioning_err_centr_rss_LR(n_idx));
fprintf('C-AOA (LR): %.4f\n', avg_positioning_err_centr_aoa_LR(n_idx));
fprintf('C-Hybrid (LR): %.4f\n', avg_positioning_err_centr_hybrid_LR(n_idx));
fprintf('C-RSS-dB (LR): %.4f\n', avg_positioning_err_centr_rss_LR_dB(n_idx));
fprintf('C-Hybrid-dB (LR): %.4f\n', avg_positioning_err_centr_hybrid_LR_dB(n_idx));

% Newly added LR variables
fprintf('D-Median (LR dB): %.4f\n', avg_positioning_err_distr_median_LR_dB(n_idx));
fprintf('C-AOA (LR CRLB): %.4f\n', avg_positioning_err_centr_aoa_LR_crlb(n_idx));
fprintf('C-Hybrid (LR CRLB): %.4f\n', avg_positioning_err_centr_hybrid_LR_crlb(n_idx));
fprintf('C-Hybrid-dB (LR CRLB): %.4f\n', avg_positioning_err_centr_hybrid_LR_dB_crlb(n_idx));
fprintf('D-Median (LR dB CRLB): %.4f\n', avg_positioning_err_distr_median_LR_dB_crlb(n_idx));

% WKNN
fprintf('C-RSS (WKNN): %.4f\n', avg_positioning_err_centr_rss_WKNN(n_idx));
fprintf('C-AOA (WKNN): %.4f\n', avg_positioning_err_centr_aoa_WKNN(n_idx));
fprintf('C-Hybrid (WKNN): %.4f\n', avg_positioning_err_centr_hybrid_WKNN(n_idx));
fprintf('C-RSS-dB (WKNN): %.4f\n', avg_positioning_err_centr_rss_WKNN_dB(n_idx));
fprintf('C-Hybrid-dB (WKNN): %.4f\n', avg_positioning_err_centr_hybrid_WKNN_dB(n_idx));

% Newly added WKNN variables
fprintf('D-Median (WKNN dB): %.4f\n', avg_positioning_err_distr_median_WKNN_dB(n_idx));
fprintf('C-AOA (WKNN CRLB): %.4f\n', avg_positioning_err_centr_aoa_WKNN_crlb(n_idx));
fprintf('C-Hybrid-dB (WKNN CRLB): %.4f\n', avg_positioning_err_centr_hybrid_WKNN_dB_crlb(n_idx));
fprintf('D-Median (WKNN dB CRLB): %.4f\n', avg_positioning_err_distr_median_WKNN_dB_crlb(n_idx));

fprintf("fcnn : %f \n",avg_positioning_err_fcnn(n_idx));

% Print average error ellipse areas
fprintf('\nAverage Error Ellipse Areas:\n');
fprintf('C-RSS: %.4f\n', avg_err_ellipse_area_centr_rss(n_idx));
fprintf('C-AOA: %.4f\n', avg_err_ellipse_area_centr_aoa(n_idx));
fprintf('C-Hybrid: %.4f\n', avg_err_ellipse_area_centr_hybrid(n_idx));
fprintf('D-Median: %.4f\n', avg_err_ellipse_area_distr_median(n_idx));
fprintf('D-Mean: %.4f\n', avg_err_ellipse_area_distr_mean(n_idx));
fprintf('D-ZS TZ 0.6: %.4f\n', avg_err_ellipse_area_distr_zscore_tz_0p6(n_idx));
fprintf('D-ZS TZ 0.8: %.4f\n', avg_err_ellipse_area_distr_zscore_tz_0p8(n_idx));
fprintf('D-ZS TZ 1: %.4f\n', avg_err_ellipse_area_distr_zscore_tz_1(n_idx));
fprintf('D-ZS TZ 1.2: %.4f\n', avg_err_ellipse_area_distr_zscore_tz_1p2(n_idx));
fprintf('D-ZS TZ 1.4: %.4f\n', avg_err_ellipse_area_distr_zscore_tz_1p4(n_idx));
fprintf('D-BAY: %.4f\n', avg_err_ellipse_area_distr_bayesian(n_idx));
fprintf('D-ZS-BAY TZ 0.6: %.4f\n', avg_err_ellipse_area_distr_zscore_tz_0p6_bayesian(n_idx));
fprintf('D-ZS-BAY TZ 0.8: %.4f\n', avg_err_ellipse_area_distr_zscore_tz_0p8_bayesian(n_idx));
fprintf('D-ZS-BAY TZ 1: %.4f\n', avg_err_ellipse_area_distr_zscore_tz_1_bayesian(n_idx));
fprintf('D-ZS-BAY TZ 1.2: %.4f\n', avg_err_ellipse_area_distr_zscore_tz_1p2_bayesian(n_idx));
fprintf('D-ZS-BAY TZ 1.4: %.4f\n', avg_err_ellipse_area_distr_zscore_tz_1p4_bayesian(n_idx));

fprintf('C-Hybrid (CRLB): %.4f\n', avg_err_ellipse_area_centr_hybrid_crlb(n_idx));
fprintf('C-AOA (CRLB): %.4f\n', avg_err_ellipse_area_centr_aoa_crlb(n_idx));
fprintf('D-Median (CRLB): %.4f\n', avg_err_ellipse_area_distr_median_crlb(n_idx));
fprintf('D-Mean (CRLB): %.4f\n', avg_err_ellipse_area_distr_mean_crlb(n_idx));
fprintf('D-ZS TZ 1 (CRLB): %.4f\n', avg_err_ellipse_area_distr_zscore_tz_1_crlb(n_idx));
fprintf('D-BAY (CRLB): %.4f\n', avg_err_ellipse_area_distr_bayesian_crlb(n_idx));
fprintf('D-ZS-BAY TZ 1 (CRLB): %.4f\n', avg_err_ellipse_area_distr_zscore_tz_1_bayesian_crlb(n_idx));

end %n_idx = 1:N_max
fprintf('\n');
fprintf('NUM_TP_INSIDE_ELLIPSE_CR: %d\n', NUM_TP_INSIDE_ELLIPSE_CR);
fprintf('NUM_TP_INSIDE_ELLIPSE_CA: %d\n', NUM_TP_INSIDE_ELLIPSE_CA);
fprintf('NUM_TP_INSIDE_ELLIPSE_CH: %d\n', NUM_TP_INSIDE_ELLIPSE_CH);
fprintf('NUM_TP_INSIDE_ELLIPSE_DD: %d\n', NUM_TP_INSIDE_ELLIPSE_DD);
fprintf('NUM_TP_INSIDE_ELLIPSE_DM: %d\n', NUM_TP_INSIDE_ELLIPSE_DM);
fprintf('NUM_TP_INSIDE_ELLIPSE_DZ: %d\n', NUM_TP_INSIDE_ELLIPSE_DZ);
fprintf('NUM_TP_INSIDE_ELLIPSE_DBAY: %d\n', NUM_TP_INSIDE_ELLIPSE_DBAY);
fprintf('NUM_TP_INSIDE_ELLIPSE_DZ_BAY: %d\n', NUM_TP_INSIDE_ELLIPSE_DZ_BAY);
fprintf('NUM_TP_INSIDE_ELLIPSE_CA_CRLB: %d\n', NUM_TP_INSIDE_ELLIPSE_CA_CRLB);
fprintf('NUM_TP_INSIDE_ELLIPSE_CH_CRLB: %d\n', NUM_TP_INSIDE_ELLIPSE_CH_CRLB);
fprintf('NUM_TP_INSIDE_ELLIPSE_DD_CRLB: %d\n', NUM_TP_INSIDE_ELLIPSE_DD_CRLB);
fprintf('NUM_TP_INSIDE_ELLIPSE_DM_CRLB: %d\n', NUM_TP_INSIDE_ELLIPSE_DM_CRLB);
fprintf('NUM_TP_INSIDE_ELLIPSE_DZ_CRLB: %d\n', NUM_TP_INSIDE_ELLIPSE_DZ_CRLB);
fprintf('NUM_TP_INSIDE_ELLIPSE_DBAY_CRLB: %d\n', NUM_TP_INSIDE_ELLIPSE_DBAY_CRLB);
fprintf('NUM_TP_INSIDE_ELLIPSE_DZ_BAY_CRLB: %d\n', NUM_TP_INSIDE_ELLIPSE_DZ_BAY_CRLB);

fprintf('\n\n');
percent_CR = NUM_TP_INSIDE_ELLIPSE_CR * 100 / (nbrOfSetups*num_tp_points*N_max);
percent_CA = NUM_TP_INSIDE_ELLIPSE_CA * 100 / (nbrOfSetups*num_tp_points*N_max);
percent_CH = NUM_TP_INSIDE_ELLIPSE_CH * 100 / (nbrOfSetups*num_tp_points*N_max);
percent_DD = NUM_TP_INSIDE_ELLIPSE_DD * 100 / (nbrOfSetups*num_tp_points*N_max);
percent_DM = NUM_TP_INSIDE_ELLIPSE_DM * 100 / (nbrOfSetups*num_tp_points*N_max);
percent_DZ = NUM_TP_INSIDE_ELLIPSE_DZ * 100 / (nbrOfSetups*num_tp_points*N_max);
percent_DBAY = NUM_TP_INSIDE_ELLIPSE_DBAY * 100 / (nbrOfSetups*num_tp_points*N_max);
percent_DZ_BAY = NUM_TP_INSIDE_ELLIPSE_DZ_BAY * 100 / (nbrOfSetups*num_tp_points*N_max);

percent_CA_CRLB = NUM_TP_INSIDE_ELLIPSE_CA_CRLB * 100 / (nbrOfSetups*num_tp_points*N_max);
percent_CH_CRLB = NUM_TP_INSIDE_ELLIPSE_CH_CRLB * 100 / (nbrOfSetups*num_tp_points*N_max);
percent_DD_CRLB = NUM_TP_INSIDE_ELLIPSE_DD_CRLB * 100 / (nbrOfSetups*num_tp_points*N_max);
percent_DM_CRLB = NUM_TP_INSIDE_ELLIPSE_DM_CRLB * 100 / (nbrOfSetups*num_tp_points*N_max);
percent_DZ_CRLB = NUM_TP_INSIDE_ELLIPSE_DZ_CRLB * 100 / (nbrOfSetups*num_tp_points*N_max);
percent_DBAY_CRLB = NUM_TP_INSIDE_ELLIPSE_DBAY_CRLB * 100 / (nbrOfSetups*num_tp_points*N_max);
percent_DZ_BAY_CRLB = NUM_TP_INSIDE_ELLIPSE_DZ_BAY_CRLB * 100 / (nbrOfSetups*num_tp_points*N_max);

% Percentages for each condition
fprintf('TPs in Ellipse CR: %.2f%%\n', percent_CR);
fprintf('TPs in Ellipse CA: %.2f%%\n', percent_CA);
fprintf('TPs in Ellipse CH: %.2f%%\n', percent_CH);
fprintf('TPs in Ellipse DD: %.2f%%\n', percent_DD);
fprintf('TPs in Ellipse DM: %.2f%%\n', percent_DM);
fprintf('TPs in Ellipse DZ: %.2f%%\n', percent_DZ);
fprintf('TPs in Ellipse DBAY: %.2f%%\n', percent_DBAY);
fprintf('TPs in Ellipse DZ_BAY: %.2f%%\n', percent_DZ_BAY);

% CRLB percentages for each condition
fprintf('TPs in Ellipse CA (CRLB): %.2f%%\n', percent_CA_CRLB);
fprintf('TPs in Ellipse CH (CRLB): %.2f%%\n', percent_CH_CRLB);
fprintf('TPs in Ellipse DD (CRLB): %.2f%%\n', percent_DD_CRLB);
fprintf('TPs in Ellipse DM (CRLB): %.2f%%\n', percent_DM_CRLB);
fprintf('TPs in Ellipse DZ (CRLB): %.2f%%\n', percent_DZ_CRLB);
fprintf('TPs in Ellipse DBAY (CRLB): %.2f%%\n', percent_DBAY_CRLB);
fprintf('TPs in Ellipse DZ_BAY (CRLB): %.2f%%\n', percent_DZ_BAY_CRLB);

diary off