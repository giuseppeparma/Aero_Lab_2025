function [centers,normals,tangents,extrema_1,extrema_2,alpha,lengths,L2G_TransfMatrix,G2L_TransfMatrix] = createPanels(geo)
% createPanels
% This function creates panels variables for Hess-Smith method
%
% INPUT
% Name                Type           Size
% geo                 struct         1x1
%
% OUTPUT  
% Name                Type           Size
% centers             matrix         Nx2
% normals             matrix         Nx2
% tangents            matrix         Nx2
% extrema_1           matrix         Nx2
% extrema_2           matrix         Nx2
% alpha               vector         Nx1
% lengths             vector         Nx1
% L2G_TransfMatrix    matrix         Nx2x2
% G2L_TransfMatrix    matrix         Nx2x2
%

    % Outputs initialitazion
    N = length(geo.x) - 1;
    
    centers = zeros(N, 2);
    normals = zeros(N, 2);
    tangents = zeros(N, 2);
    extrema_1 = zeros(N, 2);
    extrema_2 = zeros(N, 2);
    alpha = zeros(N, 1);
    L2G_TransfMatrix = zeros(N, 2, 2);
    G2L_TransfMatrix = zeros(N, 2, 2);
    
    for i = 1:N
        centers(i, 1) = (geo.x(i) + geo.x(i+1))/2;
        centers(i, 2) = (geo.y(i) + geo.y(i+1))/2;
        
        extrema_1(i, 1) = geo.x(i);
        extrema_1(i, 2) = geo.y(i);
        
        extrema_2(i, 1) = geo.x(i+1);
        extrema_2(i, 2) = geo.y(i+1);
        
        dy = extrema_2(i, 2) - extrema_1(i, 2);
        dx = extrema_2(i, 1) - extrema_1(i, 1);
        
        angle = atan2(dy, dx); % global reference system
        
        if abs(angle) < 10^(-12)
           angle = 0; 
        end
        
        alpha(i) = angle;
        sinAngle = sin(angle);
        cosAngle = cos(angle);

        if abs(sinAngle) < 10^(-12)
            sinAngle = 0;
        end
        if abs(cosAngle) < 10^(-12)
            cosAngle = 0;
        end

        L2G_TransfMatrix(i, :, :) = [cosAngle ,  -sinAngle;
                                     sinAngle,  cosAngle];
                                 
        G2L_TransfMatrix(i, :, :) = [cosAngle ,  sinAngle;
                                     -sinAngle,  cosAngle];
                                 
        normals(i, 1) = -sinAngle;
        normals(i, 2) = cosAngle;
        
        tangents(i, 1) = cosAngle;
        tangents(i, 2) = sinAngle;
        
        lengths(i) = norm(extrema_2(i, :) - extrema_1(i, :));       
    end
end