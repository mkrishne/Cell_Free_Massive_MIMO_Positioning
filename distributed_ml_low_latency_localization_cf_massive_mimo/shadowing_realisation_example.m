% Parameters
K = 5; % Number of RPs
L = 3; % Number of APs
sigma_sf = 2; % Shadowing variance
decorr = 10; % Decorrelation distance
RP_positions = linspace(0, 50, K); % RP positions (spatial locations)
TP_position = 20; % Single TP position

% Initialize correlation matrix and shadowing realizations
shadowCorrMatrix = sigma_sf^2 * ones(K, K);
shadowAPrealizations = zeros(K, L);

% Offline phase: Shadowing for RPs
for RP_idx = 1:K
    if RP_idx-1>0
        % Compute distances from the current RP to all previous RPs
        RPDistances = zeros(RP_idx-1,1);
        for i = 1:RP_idx-1
            RPDistances(i) = abs(RP_positions(RP_idx) - RP_positions(i));
        end
        
        % Compute new column for correlation matrix
        newcolumn = sigma_sf^2 * 2.^(-RPDistances / decorr);
        
        % Conditional statistics
        term1 = newcolumn' / shadowCorrMatrix(1:RP_idx-1, 1:RP_idx-1);
        meanvalues = term1 * shadowAPrealizations(1:RP_idx-1, :);
        stdvalue = sqrt(sigma_sf^2 - term1 * newcolumn);
    else
        % First RP
        %Add the UE and begin to store shadow fading correlation values
        meanvalues = 0;
        stdvalue = sigma_sf;
        newcolumn = [];
    end
    
    % Generate shadowing realizations
    rss_shadowing_offline = meanvalues + stdvalue * randn(1, L);
    shadowCorrMatrix(1:RP_idx-1,RP_idx) = newcolumn;
    shadowCorrMatrix(RP_idx,1:RP_idx-1) = newcolumn';
    shadowAPrealizations(RP_idx,:) = rss_shadowing_offline;
end

disp('Offline Shadowing Realizations (RPs):');
disp(shadowAPrealizations);

% Online phase: Shadowing for TP
TPDistance = abs(RP_positions(:) - TP_position); % Distance of TP to all RPs
newcolumn = sigma_sf^2 * 2.^(-TPDistance / decorr);
term1 = newcolumn' / shadowCorrMatrix;

tp_meanvalues = term1 * shadowAPrealizations;
tp_stdvalue = sqrt(sigma_sf^2 - term1 * newcolumn);

rss_shadowing_online = tp_meanvalues + tp_stdvalue * randn(1, L);

disp('Online Shadowing Realizations (TP):');
disp(rss_shadowing_online);
