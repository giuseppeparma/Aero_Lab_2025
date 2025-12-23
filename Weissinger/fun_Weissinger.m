function [Gamma, CL, CD_ind, alpha_induced, elapsedTimeWeissinger] = fun_Weissinger(wing, tail, U_inf_modulus, U_inf, rho)
% fun_Weissinger
% This function applies the Weissinger (optimized) method for the given
% input a/c (wing, tail and aerodynamics) data.
%
% INPUT
% Name                     Type                Size              
% wing                     struct              1x1
% tail                     struct              1x1
% U_inf_modulus            double              1x1
% U_inf                    double              1x3
% rho                      double              1x1
%
% OUTPUT
% Name                     Type                Size
% Gamma                    struct              1x1
% CL                       struct              1x1
% CD_ind                   struct              1x1
% alpha_induced            double              1xNPanelsTot
% elapsedTimeWeissinger    double              1x1
%

    tic  % starting elapsed time computation for Weissinger

    % Geometry
    components = {wing, tail};
    geometries = cell(2,1);
    for k = 1:length(components)
        comp = components{k};
        comp.controlPoints = cell(comp.ID, 1);
        comp.inducedPoints = cell(comp.ID, 1);
        comp.normals = cell(comp.ID, 1);
        comp.infiniteVortices = cell(comp.ID, 1);
        comp.vortices = cell(comp.ID, 1);
        comp.internalMesh = cell(comp.ID, 1);
        comp.extremes = cell(comp.ID, 1);
        comp.evaluationPoints = cell(comp.ID, 1);
    
        [comp.controlPoints{comp.ID}, comp.evaluationPoints{comp.ID}, comp.normals{comp.ID}, ...
         comp.infiniteVortices{comp.ID}, comp.vortices{comp.ID}, comp.internalMesh{comp.ID}] = createStructure(comp);
    
        geometries{k} = comp;
    end
    
    % Global variables definition to take into account both wing and tail
    Global.CP        = [];  % control points (3/4 chord)
    Global.EP        = [];  % evaluation point (1/4 chord)
    Global.NormalsCP = []; 
    Global.NormalsEP = [];
    Global.V1_S      = []; Global.V1_E    = [];
    Global.V_Fin_S   = []; Global.V_Fin_E = [];
    Global.V2_S      = []; Global.V2_E    = [];
    
    panelCounts = zeros(length(components), 1); 
    
    for k = 1:length(components)
        comp = geometries{k};
        N_chord = comp.chordwiseDiscr;
        N_span = comp.spanwiseSectors;
        N_panels = N_chord * N_span;
        panelCounts(k) = N_panels;
    
        tmp_CP      = zeros(N_panels, 3);
        tmp_EP      = zeros(N_span, 3);
        tmp_NormCP  = zeros(N_panels, 3);
        tmp_NormEP  = zeros(N_span, 3);
        tmp_V1S     = zeros(N_panels, 3); tmp_V1E = zeros(N_panels, 3);
        tmp_VFS     = zeros(N_panels, 3); tmp_VFE = zeros(N_panels, 3);
        tmp_V2S     = zeros(N_panels, 3); tmp_V2E = zeros(N_panels, 3);
    
        idx = 0;
        for i = 1:N_chord
            for j = 1:N_span
                idx = idx + 1;
                tmp_CP(idx,:)     = comp.controlPoints{comp.ID}{i,j}.Coords;
                tmp_NormCP(idx,:) = comp.normals{comp.ID}{i,j}.Coords;         
                tmp_V1S(idx,:)    = comp.infiniteVortices{comp.ID}{i,j}.Root.toInfty;
                tmp_V1E(idx,:)    = comp.infiniteVortices{comp.ID}{i,j}.Root.onWing;
                tmp_VFS(idx,:)    = comp.vortices{comp.ID}{i,j}.Root;
                tmp_VFE(idx,:)    = comp.vortices{comp.ID}{i,j}.Tip;
                tmp_V2S(idx,:)    = comp.infiniteVortices{comp.ID}{i,j}.Tip.onWing;
                tmp_V2E(idx,:)    = comp.infiniteVortices{comp.ID}{i,j}.Tip.toInfty;
            end
        end
    
        % Evaluation points (separated for sectors discretization)
        idx = 0;
        for j = 1:N_span
            idx = idx + 1;
            tmp_EP(idx,:)     = comp.evaluationPoints{comp.ID}{1,j}.Coords;
            tmp_NormEP(idx,:) = comp.normals{comp.ID}{1,j}.Coords; 
        end
    
        Global.CP        = [Global.CP;        tmp_CP];
        Global.EP        = [Global.EP;        tmp_EP];
        Global.NormalsCP = [Global.NormalsCP; tmp_NormCP];
        Global.NormalsEP = [Global.NormalsEP; tmp_NormEP];
        Global.V1_S      = [Global.V1_S;      tmp_V1S];
        Global.V1_E      = [Global.V1_E;      tmp_V1E];
        Global.V_Fin_S   = [Global.V_Fin_S;   tmp_VFS];
        Global.V_Fin_E   = [Global.V_Fin_E;   tmp_VFE];
        Global.V2_S      = [Global.V2_S;      tmp_V2S];
        Global.V2_E      = [Global.V2_E;      tmp_V2E];
    end
    
    N_PanelsTotal = sum(panelCounts);
    
    % Linear system: complete for two surfaces
    % A matrix
    A = zeros(N_PanelsTotal, N_PanelsTotal);
    for i = 1:N_PanelsTotal
        ControlPoint = Global.CP(i, :);
        Normal = Global.NormalsCP(i, :)';
    
        % Velocity per unit circulation
        U_inf1 = vortexInfluence(ControlPoint, Global.V1_S,    Global.V1_E);
        U_fin  = vortexInfluence(ControlPoint, Global.V_Fin_S, Global.V_Fin_E);
        U_inf2 = vortexInfluence(ControlPoint, Global.V2_S,    Global.V2_E);
    
        U_tot = U_inf1 + U_fin + U_inf2;
    
        % A matrix
        A(i, :) = (U_tot * Normal)';
    end
    
    % b known term
    b = -dot(repmat(U_inf, N_PanelsTotal, 1), Global.NormalsCP, 2);
    
    solution = A \ b;  % solution of the linear system (circulation)
    
    % Circulation
    idx_Start_Wing = 1;                idx_End_Wing = panelCounts(1);
    idx_Start_Tail = panelCounts(1)+1; idx_End_Tail = N_PanelsTotal;
    
    Gamma_Wing_Vec = solution(idx_Start_Wing:idx_End_Wing);
    Gamma_Tail_Vec = solution(idx_Start_Tail:idx_End_Tail);
    
    % Reshape in grid
    Wing_Grid_Gamma = reshape(Gamma_Wing_Vec, [geometries{1}.spanwiseSectors, geometries{1}.chordwiseDiscr])';
    Tail_Grid_Gamma = reshape(Gamma_Tail_Vec, [geometries{2}.spanwiseSectors, geometries{2}.chordwiseDiscr])';

    Gamma.Wing_Vec = Gamma_Wing_Vec;
    Gamma.Tail_Vec = Gamma_Tail_Vec;
    Gamma.Wing_Grid_Gamma = Wing_Grid_Gamma;
    Gamma.Tail_Grid_Gamma = Tail_Grid_Gamma;

    % figure('Name','Circulation Distribution',NumberTitle='off');
    % subplot(2, 1, 1); 
    % pcolor(Tail_Grid_Gamma);
    % title('Tail Circulation');
    % shading interp; 
    % colorbar;
    % 
    % subplot(2, 1, 2); 
    % pcolor(Wing_Grid_Gamma); 
    % title('Wing Circulation');
    % shading interp; 
    % colorbar; 
    
    % Lift 2D and 3D of the wing
    % Kutta-Joukowsky theorem    
    Gamma_Span_Dist_wing = sum(Wing_Grid_Gamma, 1);
    Gamma_Span_Dist_tail = sum(Tail_Grid_Gamma, 1);
    
    L_2D_wing = rho * U_inf_modulus .* Gamma_Span_Dist_wing .* cosd(wing.dihedralAngle);
    L_2D_tail = rho * U_inf_modulus .* Gamma_Span_Dist_tail .* cosd(tail.dihedralAngle);
    L_2D = [L_2D_wing, L_2D_tail];
    
    CL_2D_wing = L_2D_wing ./ (0.5 * rho * U_inf_modulus^2 * wing.MAC);
    CL_2D_tail = L_2D_tail ./ (0.5 * rho * U_inf_modulus^2 * tail.MAC);
    
    L_3D_wing = wing.delta_b * sum(L_2D_wing);
    L_3D_tail = tail.delta_b * sum(L_2D_tail);
    L_3D = L_3D_wing + L_3D_tail;
    
    CL_3D = L_3D / (0.5 * rho * U_inf_modulus^2 * wing.surface);

    CL.wing2D = CL_2D_wing;
    CL.tail2D = CL_2D_tail;
    CL.tot3D  = CL_3D;
    
    % Induced Drag 2D and 3D of the wing
    % Computed from induced velocity and angle
    for k = 1:length(components)
        comp = geometries{k};
        N_span = comp.spanwiseSectors;
        panelCounts(k) = N_span;
    end 
    
    alpha_induced = zeros(1, wing.spanwiseSectors);
    for i = 1:sum(panelCounts)
        EvalPoint = Global.EP(i, :);
        Normal    = Global.NormalsEP(i, :)';
    
        % Semi-infinite vortices contribution
        U_ind_inf1 = vortexInfluence(EvalPoint, Global.V1_S, Global.V1_E);
        U_ind_inf2 = vortexInfluence(EvalPoint, Global.V2_S, Global.V2_E);
    
        U_ind_tot = sum((U_ind_inf1 + U_ind_inf2).*solution, 1);
        alpha_induced(1,i) = atan( dot(U_ind_tot, Normal) / U_inf_modulus );
    end
    
    D_ind_2D = L_2D .* sin(alpha_induced);
    D_ind_2D_wing = D_ind_2D(1:wing.spanwiseSectors);
    D_ind_2D_tail = D_ind_2D((wing.spanwiseSectors + 1) : (wing.spanwiseSectors + tail.spanwiseSectors));
    
    D_ind_3D_wing = wing.delta_b * sum(abs(D_ind_2D_wing));
    D_ind_3D_tail = tail.delta_b * sum(abs(D_ind_2D_tail));
    D_ind_3D = D_ind_3D_wing + D_ind_3D_tail;
    
    CD_ind_3D = D_ind_3D / (0.5 * rho * U_inf_modulus^2 * wing.surface);

    CD_ind.tot3D = CD_ind_3D;

    % Elapsed time for the implemented Weissinger method (not including plots)
    elapsedTimeWeissinger = toc;
end