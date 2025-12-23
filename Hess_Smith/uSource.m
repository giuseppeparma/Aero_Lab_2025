function U = uSource(center, extreme_1, extreme_2, L2G_TransfMatrix, G2L_TransfMatrix)
% uSource
% This function computes the source induced velocity
%
% INPUT 
% Name                Type           Size
% center              scalar         1x1
% extreme_1           scalar         1x1
% extreme_2           scalar         1x1
% L2G_TransfMatrix    matrix         2x2
% G2L_TransfMatrix    matrix         2x2
%
% OUTPUT
% Name                Type           Size
% U                   vector         2x1
%

    % Global to local coordinates   
    center = G2L_TransfMatrix * center;
    extreme_1 = G2L_TransfMatrix * extreme_1;
    extreme_2 = G2L_TransfMatrix * extreme_2;
    
    % Local u and v
    r1 = center - extreme_1;
    theta_1 = atan2(r1(2), r1(1));
    r2 = center - extreme_2;
    theta_2 = atan2(r2(2), r2(1));
    
    if abs(theta_1) < 10^(-12) && abs(theta_2) > 3
        theta_1 = 0; 
        theta_2 = pi; 
    elseif abs(theta_2) < 10^(-12) && abs(theta_1) > 3
        theta_2 = 0;
        theta_1 = -pi; 
    end

    u = -(0.5/pi) * log(norm(r2)/norm(r1));
    v = theta_2 - theta_1;
    v = v / (2*pi);
    
    % Local to global coordinates
    U = L2G_TransfMatrix * [u;v];
    if abs(U(1)) < 10^(-12)
        U(1) = 0;
    end
    if abs(U(2)) < 10^(-12)
        U(2) = 0;
    end
end