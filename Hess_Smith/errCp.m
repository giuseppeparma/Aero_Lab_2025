function Cp_err = errCp(NPanels,centers,x_Cp,Cp_HS,Cp_XF)
% errCp
% This function computes the error for the Cp between Hess-Smith and XFoil
%
% INPUT
% Name           Type           Size
% NPanels        scalar         1x1
% centers        matrix         Nx2
% x_Cp           vector         (N+1)x1
% Cp_HS          vector         Nx1
% Cp_XF          vector         (N+1)x1
%
% OUTPUT
% Name           Type           Size
% Cp_err         vector         (N+1)x1
%
    % Interpolation
    p = polyfit(centers(:,1),Cp_HS,NPanels+1);
    Cp_HS = polyval(p,x_Cp);

    Cp_err = abs(Cp_HS - Cp_XF);
end