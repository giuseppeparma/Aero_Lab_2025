function U = vortexInfluence(targetPoint, startPoint, endPoint)
% vortexInfluence
% This function creates the unit-circulation velocity due to vortices.
%
% INPUT
% Name               Type                Size              
% targetPoint        double              1x3
% startPoint         double              1x3
% endPoint           double              1x3
%
% OUTPUT
% Name               Type                Size
% U                  double              1x3
% 

    tol = 1e-8;

    r1 = targetPoint - startPoint;
    r2 = targetPoint - endPoint;
    r0 = endPoint - startPoint;

    r1_modulus = sqrt( sum(r1.^2, 2) );
    r2_modulus = sqrt( sum(r2.^2, 2) );

    r1_x_r2 = cross(r1, r2, 2);
    r1_x_r2_sq = sum(r1_x_r2.^2, 2);

    if r1_x_r2_sq < tol
       r1_x_r2_sq = tol;
    end

    dot_term = (sum(r0 .* r1, 2) ./ r1_modulus) - (sum(r0 .* r2, 2) ./ r2_modulus);
    scalar_part = (1/(4*pi) .* dot_term) ./ r1_x_r2_sq;

    U = r1_x_r2 .* scalar_part;
end
