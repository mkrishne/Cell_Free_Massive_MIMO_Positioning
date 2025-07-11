function [RP_positions,TP_positions,AP_positions] = cell_free_layout_setup_random_RP(RP_positions_per_row,squareLength,L,minDistanceUE2AP,minDistanceAP2AP,num_tp_points)
K = RP_positions_per_row^2;

RP_positions_X = rand(1, K) * squareLength;
RP_positions_Y = rand(1, K) * squareLength;
RP_positions = RP_positions_X + 1i*RP_positions_Y;

%initialize infinite distance to ensure second condition inside for loop is satisfied at start
AP_positions = 2*squareLength*(ones(1,L) + 1i*ones(1,L)); 
%Random AP locations with uniform distribution


for AP_idx = 1:L
    AP_pos_rand = (rand + 1i*rand) * squareLength;
    %if the dist is less than 10m between any two APs or RP positions and any AP then reposition the AP
    while( (min(abs(RP_positions - AP_pos_rand), [], "all") < minDistanceUE2AP) || ...
           (min(abs(AP_positions - AP_pos_rand), [], "all") < minDistanceAP2AP) ) 
        AP_pos_rand = (rand + 1i*rand) * squareLength;
    end
    AP_positions(AP_idx) =  AP_pos_rand;
end

TP_positions = zeros(1,num_tp_points);
for TP_idx = 1:num_tp_points
    TP_pos_rand = (rand + 1i*rand) * squareLength;
    while(min(abs(AP_positions - TP_pos_rand), [], "all") < minDistanceUE2AP)
        TP_pos_rand = (rand + 1i*rand) * squareLength;
    end
    TP_positions(TP_idx) = TP_pos_rand;
end




