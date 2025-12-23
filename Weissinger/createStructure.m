function [controlPoints, evaluationPoints, normals, infiniteVortices, vortices, internalMesh] = createStructure(config)
% createStructure
% This function creates the structure for the input lifting surface.
%
% INPUT
% Name                  Type                Size              
% config                struct              1x1
%
% OUTPUT
% Name                  Type                Size
% controlPoints         cell                1x1
% evaluationPoints      cell                1x1
% normals               cell                1x1
% infiniteVortices      cell                1x1
% vortices              cell                1x1
% internalMesh          cell                1x1
%

    N_semispan = config.semispanwiseDiscr;
    N_span = 2 * N_semispan;
    N_chord = config.chordwiseDiscr;   
    
    YAngle = config.rotationangleY;
    YRot = [ cosd(YAngle)   0       sind(YAngle); 
                0           1           0; 
            -sind(YAngle)   0       cosd(YAngle)];

    XAngle = config.dihedralAngle;
    XRot = [    1           0           0; 
                0       cosd(XAngle) -sind(XAngle); 
                0       sind(XAngle)  cosd(XAngle)];
    
    
    % Half-wing
    RootChord = config.rootChord;
    TipChord  = config.tipChord;
    SemiSpan  = config.semispan;
    Sweep     = config.sweepAngle;

    P_RootLE = [0, 0, 0];
    P_RootTE = (YRot * [RootChord; 0; 0])'; 
    
    % Adding sweep and dihedral angles
    dx_sweep = RootChord/4 + SemiSpan * tand(Sweep) - TipChord/4;
    P_TipLE_base = [dx_sweep, SemiSpan, 0];

    % Adding rotations
    P_TipLE = (YRot * (XRot * P_TipLE_base'))';
    P_TipTE = P_TipLE + (YRot * [TipChord; 0; 0])';

    % Positioning
    Offset = [config.LEpositionX, config.LEpositionY, config.LEpositionZ];
    P_RootLE = P_RootLE + Offset; P_RootTE = P_RootTE + Offset;
    P_TipLE  = P_TipLE + Offset;  P_TipTE  = P_TipTE + Offset;

    % Mesh
    eta = linspace(0, 1, N_semispan + 1);  % spanwise
    xi  = linspace(0, 1, N_chord + 1)';      % chordwise (column linspace)

    LE_vec = P_RootLE + (P_TipLE - P_RootLE) .* eta'; 
    TE_vec = P_RootTE + (P_TipTE - P_RootTE) .* eta';
   
    LE_grid = repmat(reshape(LE_vec, [1, N_semispan+1, 3]), [N_chord+1, 1, 1]);
    TE_grid = repmat(reshape(TE_vec, [1, N_semispan+1, 3]), [N_chord+1, 1, 1]);
    Xi_grid = repmat(xi, [1, N_semispan+1, 3]);

    Mesh_Right = LE_grid + Xi_grid .* (TE_grid - LE_grid);
    
    % Mirroring for the other half-wing
    Mesh_Left = Mesh_Right;
    Mesh_Left(:,:,2) = -Mesh_Left(:,:,2); 

    Mesh_Left = flip(Mesh_Left, 2);

    % Complete mesh
    Mesh_Total = cat(2, Mesh_Left(:, 1:end-1, :), Mesh_Right);
    
    % Values for root and tip LE and TE
    P1 = Mesh_Total(1:end-1, 1:end-1, :); % root LE
    P2 = Mesh_Total(1:end-1, 2:end, :);   % tip LE
    P3 = Mesh_Total(2:end, 2:end, :);     % tip TE
    P4 = Mesh_Total(2:end, 1:end-1, :);   % root TE
    
    % Quarter chord position root and tip (for finite vortices)
    Q_Root = P1 + 0.25*(P4 - P1);
    Q_Tip  = P2 + 0.25*(P3 - P2);

    % Three-quarter chord position root and tip (for control points)\
    TQ_Root = P1 + 0.75*(P4 - P1);
    TQ_Tip  = P2 + 0.75*(P3 - P2);
    
    Vortex_Roots = Q_Root;
    Vortex_Tips  = Q_Tip;
    CPs = (TQ_Root + TQ_Tip) / 2;
    
    % Normals
    VecA = P3 - P4; % TE Tip - TE Root (spanwise at TE)
    VecB = P1 - P4; % LE Root - TE Root (chordwise between LE and TE)
    NormalsVec = cross(VecA, VecB, 3);
    NormalsModulus = sqrt(sum(NormalsVec.^2, 3));
    NormalsVec = NormalsVec ./ NormalsModulus;
    
    % Semi-infinite definitions
    InfLength = 50 * config.rootChord;
    InfDist = [InfLength, 0, 0];
    InfDistGrid = repmat(reshape(InfDist, [1,1,3]), [N_chord, N_span, 1]);
    
    % Output development
    controlPoints    = cell(N_chord, N_span);
    normals          = cell(N_chord, N_span);
    infiniteVortices = cell(N_chord, N_span);
    vortices         = cell(N_chord, N_span);
    internalMesh     = cell(N_chord, N_span);
    
    for i = 1:N_chord
        for j = 1:N_span
            controlPoints{i,j}.Coords = squeeze(CPs(i,j,:))';
            normals{i,j}.Coords = squeeze(NormalsVec(i,j,:))';
            vortices{i,j}.Root = squeeze(Vortex_Roots(i,j,:))';
            vortices{i,j}.Tip  = squeeze(Vortex_Tips(i,j,:))';
            infiniteVortices{i,j}.Root.onWing  = squeeze(Vortex_Roots(i,j,:))';
            infiniteVortices{i,j}.Root.toInfty = squeeze(Vortex_Roots(i,j,:)+InfDistGrid(i,j,:))';
            infiniteVortices{i,j}.Tip.onWing   = squeeze(Vortex_Tips(i,j,:))';
            infiniteVortices{i,j}.Tip.toInfty  = squeeze(Vortex_Tips(i,j,:)+InfDistGrid(i,j,:))';
            
            internalMesh{i,j}.LERoot = squeeze(P1(i,j,:))';
            internalMesh{i,j}.LEtip  = squeeze(P2(i,j,:))';
            internalMesh{i,j}.TEtip  = squeeze(P3(i,j,:))';
            internalMesh{i,j}.TERoot = squeeze(P4(i,j,:))';
        end
    end

    % Evaluation points (separated for sectors discretization)
    evaluationPoints = cell(1, N_span);
    for j = 1:N_span
        % Geometric mean LE and TE of the j-th sector
        SecLE = (internalMesh{1,j}.LERoot + internalMesh{1,j}.LEtip)/2;
        SecTE = (internalMesh{end,j}.TERoot + internalMesh{end,j}.TEtip)/2;
        
        Ep = zeros(1,3);
        Ep(1) = SecLE(1) + (SecTE(1) - SecLE(1))/4;
        Ep(2) = SecLE(2);
        Ep(3) = SecLE(3) + (SecTE(3) - SecLE(3))/4;
        evaluationPoints{1,j}.Coords = Ep;
    end
end