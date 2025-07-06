function crlb = find_crlb(Nsnapshots,array_spacing,ang_in,half_angular_spread,steeringvec,signal_pow_per_ant,R)

    num_ant_elem = size(steeringvec,1);
    eta =2*pi*array_spacing*deg2rad(half_angular_spread)*cosd(ang_in);
    elem_spacing_vec = (0:num_ant_elem-1).' - (0:num_ant_elem-1);
    J_0 = besselj(0,eta*elem_spacing_vec);
    J_1 = besselj(1,eta*elem_spacing_vec);
    J_2 = besselj(2,eta*elem_spacing_vec);
    J_3 = besselj(3,eta*elem_spacing_vec);

    steeringvec_matrix = signal_pow_per_ant*(steeringvec*steeringvec');
    del_R_term1 = (-1i*2*pi*array_spacing*cosd(ang_in))*((J_0 + J_2).* elem_spacing_vec .* steeringvec_matrix);
    del_R_term2 = (2*pi*array_spacing*deg2rad(half_angular_spread)*sind(ang_in))*((0.5*J_1 + 0.5*J_3).* elem_spacing_vec .* steeringvec_matrix);

    %del_R_term1 = (1i*2*pi*array_spacing*cosd(ang_in))*(J_0 .* elem_spacing_vec .* steeringvec_matrix);
    %del_R_term2 = (2*pi*array_spacing*deg2rad(half_angular_spread)*sind(ang_in))*(J_1 .* elem_spacing_vec .* steeringvec_matrix);
    del_R = del_R_term1 + del_R_term2;
    inv_R = inv(R);
    F = trace(inv_R*del_R*inv_R*del_R);
    F = real(F);
    crlb = 1/(F*Nsnapshots);
end
