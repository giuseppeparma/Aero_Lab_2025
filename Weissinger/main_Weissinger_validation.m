clear
close all
clc

tic

%% ---- Weissinger method ----
% Main used for validation: simple rectangular wing is used to compute the
%                           aerodynamic parameters for validation purpose.

% Data
U_inf_modulus = 1;  % [m/s]
beta = 0;  % [deg]
U_inf = [cosd(beta) sind(beta) 0] .* U_inf_modulus;
rho = 1.225;  % [Kg/m^3]

% Wing structure data
wing.ID = 1;

wing.dihedralAngle = 0;  % [deg]
wing.sweepAngle = 0;  % [deg]

wing.aspectRatio = 8;        
wing.taperRatio = 1;         
wing.span = 8;               
wing.semispan = wing.span/2;
wing.rootChord = 1;
wing.tipChord = wing.rootChord*wing.taperRatio;
wing.MAC = 2/3*wing.rootChord*(1+wing.taperRatio+wing.taperRatio^2)/(1+wing.taperRatio);
wing.surface = 2*(wing.semispan*wing.rootChord*(1+wing.taperRatio)/2);
wing.surfaceProjected = wing.surface*cosd(wing.dihedralAngle);

wing.semispanwiseDiscr = 20;
wing.chordwiseDiscr = 20;

wing.spanwiseSectors = 2*wing.semispanwiseDiscr;

wing.LEpositionX = 0; 
wing.LEpositionY = 0; 
wing.LEpositionZ = 0;

wing.rotationangleX = 0;  % [deg]
wing.rotationangleY = 2;  % [deg] (AoA)
wing.rotationangleZ = 0;  % [deg]

% Geometry
wing.controlPoints = cell(wing.ID, 1);
wing.normals = cell(wing.ID, 1);
wing.infiniteVortices = cell(wing.ID, 1);
wing.vortices = cell(wing.ID, 1);
wing.internalMesh = cell(wing.ID, 1);

[wing.controlPoints{wing.ID}, wing.evaluationPoints{wing.ID}, wing.normals{wing.ID}, ...
 wing.infiniteVortices{wing.ID}, wing.vortices{wing.ID}, wing.internalMesh{wing.ID}] = createStructure(wing);

N_chord = wing.chordwiseDiscr;  % chordwise panels number
N_span  = 2*wing.semispanwiseDiscr; % spanwise panels number
NPanelsTot = N_chord * N_span;

CP_coords      = zeros(NPanelsTot, 3);  % control points
CP_normals = zeros(NPanelsTot, 3);  % normals

V_Inf_1_Start  = zeros(NPanelsTot, 3); V_Inf_1_End = zeros(NPanelsTot, 3);  % first semi-infinite vortices (root)
V_Fin_Start    = zeros(NPanelsTot, 3); V_Fin_End   = zeros(NPanelsTot, 3);  % finite vortices
V_Inf_2_Start  = zeros(NPanelsTot, 3); V_Inf_2_End = zeros(NPanelsTot, 3);  % second semi-infinite vortices (tip)

idx = 0;
for i = 1:N_chord
    for j = 1:N_span
        idx = idx + 1;
        
        CP_coords(idx, :)      = wing.controlPoints{wing.ID}{i, j}.Coords;
        CP_normals(idx, :)     = wing.normals{wing.ID}{i, j}.Coords;
        
        % First semi-infinite vortex (root)
        V_Inf_1_Start(idx, :)  = wing.infiniteVortices{wing.ID}{i, j}.Root.toInfty;
        V_Inf_1_End(idx, :)    = wing.infiniteVortices{wing.ID}{i, j}.Root.onWing;
        
        % Finite vortex
        V_Fin_Start(idx, :)    = wing.vortices{wing.ID}{i, j}.Root;
        V_Fin_End(idx, :)      = wing.vortices{wing.ID}{i, j}.Tip;
        
        % Second semi-infinite vortex (tip)
        V_Inf_2_Start(idx, :)  = wing.infiniteVortices{wing.ID}{i, j}.Tip.onWing;
        V_Inf_2_End(idx, :)    = wing.infiniteVortices{wing.ID}{i, j}.Tip.toInfty;
    end
end

% Linear system
% A matrix
A = zeros(NPanelsTot, NPanelsTot);  % A matrix
for i = 1:NPanelsTot
    ControlPoint = CP_coords(i, :);
    Normal       = CP_normals(i, :)'; 
    
    % Velocity per unit circulation
    U_inf1  = vortexInfluence(ControlPoint, V_Inf_1_Start, V_Inf_1_End);
    U_fin   = vortexInfluence(ControlPoint, V_Fin_Start,   V_Fin_End);
    U_inf2  = vortexInfluence(ControlPoint, V_Inf_2_Start, V_Inf_2_End);
   
    U_tot   = U_inf1 + U_fin + U_inf2;
   
    % A matrix
    A(i, :) = (U_tot * Normal)'; 
end

% b known term
b = -dot(repmat(U_inf, NPanelsTot, 1), CP_normals, 2);

solution = A \ b;  % solution of the linear system (circulation)

wing.circulation{wing.ID} = zeros(N_chord, N_span);
idx = 0;
for i = 1:N_chord
    for j = 1:N_span
        idx = idx + 1;
        wing.circulation{wing.ID}(i, j) = solution(idx);
    end
end

figure('Name','Circulation distribution',NumberTitle='off');
pcolor(wing.circulation{wing.ID});
title('Circulation on wing panels (\Gamma_{w})');
xlabel('Spanwise'); 
ylabel('Chordwise');
shading interp
colorbar;

% Lift 2D and 3D of the wing
% Kutta-Joukowsky theorem
delta_b = wing.span / wing.spanwiseSectors; 
wing.circulationSector = sum(wing.circulation{wing.ID}, 1)'; 

L_2D = rho * U_inf_modulus * wing.circulationSector .* cosd(wing.dihedralAngle);
CL_2D = L_2D ./ (0.5 * rho * U_inf_modulus^2 * wing.MAC);

L_3D = delta_b * sum(L_2D);
CL_3D = L_3D / (0.5 * rho * U_inf_modulus^2 * wing.surface);

span_centers = (-wing.span/2 + delta_b/2):delta_b:(wing.span/2 - delta_b/2);
span_full = linspace(-wing.span/2, wing.span/2, wing.spanwiseSectors);

figure('Name','Lift coefficient distribution',NumberTitle='off');
subplot(2, 1, 1);
plot(span_centers, CL_2D, 'or');
title(['C_{l} at sector centers (Total C_{L} = ' num2str(CL_3D, '%.4f') ')']);
xlabel('span [m]'); 
ylabel('C_{l}'); 
grid on;

subplot(2, 1, 2);
plot(span_full, CL_2D, 'b', 'LineWidth', 1.5);
title('Spanwise distribution interpolated');
xlabel('span [m]'); 
ylabel('C_{l}'); 
grid on;

% Induced Drag 2D and 3D of the wing
% Computed from induced velocity and angle
EP_coords = zeros(wing.spanwiseSectors, 3);  % evaluation points
EP_normals = zeros(wing.spanwiseSectors, 3);  % normals (for sectors)

for k = 1:wing.spanwiseSectors
    EP_coords(k, :) = wing.evaluationPoints{wing.ID}{1, k}.Coords;
    EP_normals(k, :) = wing.normals{wing.ID}{1, k}.Coords; 
end

U_induced = zeros(wing.spanwiseSectors, 3);
alpha_induced = zeros(wing.spanwiseSectors, 1);
wing.circulation = solution; 
for i = 1:wing.spanwiseSectors
    EvalPoint = EP_coords(i, :);
    e_i = EP_normals(i, :);
    
    % Semi-infinite vortices contribution
    u1 = vortexInfluence(EvalPoint, V_Inf_1_Start, V_Inf_1_End);
    u2 = vortexInfluence(EvalPoint, V_Inf_2_Start, V_Inf_2_End);
    
    U_ind_contrib = (u1 + u2) .* wing.circulation;  
    U_ind_total = sum(U_ind_contrib, 1);
    
    U_induced(i, :) = U_ind_total;
    alpha_induced(i) = atan( dot(U_ind_total, e_i) / U_inf_modulus );
end

D_ind_2D = L_2D .* abs(sin(alpha_induced));
D_ind_3D = delta_b * sum(D_ind_2D);
CD_ind_3D = D_ind_3D / (0.5 * rho * U_inf_modulus^2 * wing.surface);

figure('Name','Induced angle distribution',NumberTitle='off');
subplot(2, 1, 1);
plot(span_centers, rad2deg(alpha_induced), 'or', 'MarkerFaceColor', 'r');
title(['Induced angle (Total C_{Di} = ' num2str(CD_ind_3D, '%.5f') ')']);
xlabel('Span y [m]'); 
ylabel('\alpha_{ind} [deg]'); 
grid on;

subplot(2, 1, 2);
plot(span_full, rad2deg(alpha_induced), 'b', 'LineWidth', 1.5);
title('Induced angle distribution interpolated');
xlabel('Span y [m]'); 
ylabel('\alpha_{ind} [deg]'); 
grid on;

% Elapsed time for computation
elapsedTime = toc;
fprintf('Elapsed time as consequence of optimization: %.4f [s]\n\n', elapsedTime);