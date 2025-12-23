function Cp = pressureCoeff (NPanels,centers,extrema_1,extrema_2,tangents,L2G_TransfMatrix,G2L_TransfMatrix,solution,U_inf)
% pressureCoeff
% This function computes the pressure coefficient
%
% INPUT 
% Name                Type           Size
% Npanels             scalar         1x1
% centers             matrix         Nx2
% extrema_1           matrix         Nx2
% extrema_2           matrix         Nx2
% tangents            matrix         Nx2
% L2G_TransfMatrix    matrix         Nx2x2
% G2L_TransfMatrix    matrix         Nx2x2
% solution            vector         (N+1)x1
% U_inf               vector         2x1
%
% OUTPUT
% Name                Type           Size
% Cp                  vector         Nx1
%

    U_tot = zeros(NPanels,1);
    V_tot = zeros(NPanels,1);
    Cp = zeros(NPanels,1);

    for i = 1:NPanels
         source_induced_tot = [0; 0];
         vortex_induced_tot = [0; 0];
             for j = 1:NPanels
                Temp_L2G_TransfMatrix = squeeze(L2G_TransfMatrix(j, :, :));
                Temp_G2L_TransfMatrix = squeeze(G2L_TransfMatrix(j, :, :));

                source_panel = solution(j)*uSource(centers(i,:)',extrema_1(j,:)',extrema_2(j,:)',Temp_L2G_TransfMatrix,Temp_G2L_TransfMatrix);
                source_induced_tot = source_induced_tot + source_panel;

                vortex_panel = solution(NPanels + 1)*uVortex(centers(i,:)',extrema_1(j,:)',extrema_2(j,:)',Temp_L2G_TransfMatrix,Temp_G2L_TransfMatrix);
                vortex_induced_tot = vortex_induced_tot + vortex_panel;
             end

       U_tot(i) = U_inf(1) + source_induced_tot(1) + vortex_induced_tot(1);
       V_tot(i) = U_inf(2) + source_induced_tot(2) + vortex_induced_tot(2);
             
       U_pert = [U_tot(i); V_tot(i)];
       Cp(i) = 1 - (dot(U_pert,tangents(i,:))^2)/(norm(U_inf))^2;
    end
end