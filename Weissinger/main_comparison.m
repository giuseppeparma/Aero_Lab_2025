clear
close all
clc

%% ---- Cessna 172 SkyHawk ----
% Weissinger method is applied by calling its function after defining the
% preliminary aerodynamics and wing and tail structures parameters.

% Aerodynamic data
U_inf_modulus = 20; % [m/s]
beta = 0;  % [deg]
U_inf = [cosd(beta) sind(beta) 0] .* U_inf_modulus;
rho = 1.225;  % [Kg/m^3]

% Wing structure data
Cessna.wing.ID = 1;

Cessna.wing.dihedralAngle = 1.73;  % [deg]   
Cessna.wing.sweepAngle = 1.32;  % [deg] 

Cessna.wing.aspectRatio = 7.4895;
Cessna.wing.taperRatio = 0.6875;
Cessna.wing.span = 11.0032;
Cessna.wing.semispan = Cessna.wing.span/2;
Cessna.wing.rootChord = 1.6256;
Cessna.wing.tipChord = Cessna.wing.rootChord * Cessna.wing.taperRatio;
Cessna.wing.MAC = 2/3*Cessna.wing.rootChord*(1 + Cessna.wing.taperRatio + Cessna.wing.taperRatio^2)/(1 + Cessna.wing.taperRatio);
Cessna.wing.surface = 2*(Cessna.wing.semispan*Cessna.wing.rootChord*(1 + Cessna.wing.taperRatio)/2);
Cessna.wing.surfaceProjected = Cessna.wing.surface*cosd(Cessna.wing.dihedralAngle);

Cessna.wing.semispanwiseDiscr = 20; 
Cessna.wing.chordwiseDiscr = 20;

Cessna.wing.spanwiseSectors = 2 * Cessna.wing.semispanwiseDiscr;
Cessna.wing.delta_b = Cessna.wing.span / Cessna.wing.spanwiseSectors;
   
Cessna.wing.LEpositionX = 0; 
Cessna.wing.LEpositionY = 0; 
Cessna.wing.LEpositionZ = 0;

Cessna.wing.rotationangleX = 0;  % [deg]
% Cessna.wing.rotationangleY = alpha;  defined as sequence for the Weissinger function
Cessna.wing.rotationangleZ = 0;  % [deg]

% Tail structure data
Cessna.tail.ID = 2;

Cessna.tail.dihedralAngle = 0;  % [deg]
Cessna.tail.sweepAngle = 0;  % [deg]

Cessna.tail.aspectRatio = 4.5; 
Cessna.tail.taperRatio = 0.6667;
Cessna.tail.span = 3.45; 
Cessna.tail.semispan = Cessna.tail.span/2;
Cessna.tail.rootChord = 1.05;
Cessna.tail.tipChord = Cessna.tail.rootChord * Cessna.tail.taperRatio;
Cessna.tail.MAC = 2/3*Cessna.tail.rootChord*(1 + Cessna.tail.taperRatio + Cessna.tail.taperRatio^2)/(1 + Cessna.tail.taperRatio);
Cessna.tail.surface = 2*(Cessna.tail.semispan*Cessna.tail.rootChord*(1 + Cessna.tail.taperRatio)/2);
Cessna.tail.surfaceProjected = Cessna.tail.surface*cosd(Cessna.tail.dihedralAngle);

Cessna.tail.semispanwiseDiscr = 10;
Cessna.tail.chordwiseDiscr = 10;

Cessna.tail.spanwiseSectors = 2 * Cessna.tail.semispanwiseDiscr;
Cessna.tail.delta_b = Cessna.tail.span / Cessna.tail.spanwiseSectors;

Cessna.tail.LEpositionX = 4.5;  % Positive to have the tail behind the wing
Cessna.tail.LEpositionY = 0; 
Cessna.tail.LEpositionZ = -0.5;  % Negative to have the tail below the wing

Cessna.tail.rotationangleX = 0;  % [deg]
Cessna.tail.rotationangleY = -2;  % [deg] (tail incidence)
Cessna.tail.rotationangleZ = 0;  % [deg]

alpha_seq = -4:2:15;  % [deg]
Cessna.CL_3D = zeros(length(alpha_seq), 1);
Cessna.CD_ind_3D = zeros(length(alpha_seq), 1);
for i = 1:length(alpha_seq)
    Cessna.wing.rotationangleY = alpha_seq(i);
    [Gamma, CL, CD_ind, alpha_induced, elapsedTimeWeissinger] = fun_Weissinger(Cessna.wing, Cessna.tail, U_inf_modulus, U_inf, rho);
    
    % Data for polar plot (for AoA given range)
    Cessna.CL_3D(i, 1) = CL.tot3D;
    Cessna.CD_ind_3D(i, 1) = CD_ind.tot3D;
    
    % Data for circulation, CL distribution and induced alpha plot (at specific AoA)
    if Cessna.wing.rotationangleY <= 3
       Cessna.Wing_Grid_Gamma = Gamma.Wing_Grid_Gamma;
       Cessna.Tail_Grid_Gamma = Gamma.Tail_Grid_Gamma;

       Cessna.CL_2D_wing = CL.wing2D;
       Cessna.CL_2D_tail = CL.tail2D;

       Cessna.alpha_induced = alpha_induced;
    end
    
    % Elapsed time at first Weissinger function calls
    if Cessna.wing.rotationangleY == alpha_seq(1)
       elapsedTimeWeissinger_first = elapsedTimeWeissinger;
    end
end


%% ---- Tecnam P2002-JF ----
% Weissinger method is applied by calling its function after defining the
% preliminary aerodynamics and wing and tail structures parameters.

% Wing structure data
Tecnam.wing.ID = 1;

Tecnam.wing.dihedralAngle = 5.00;  % [deg]   
Tecnam.wing.sweepAngle = 0;  % [deg] 

Tecnam.wing.aspectRatio = 6.43;
Tecnam.wing.taperRatio = 0.66;
Tecnam.wing.span = 8.60;
Tecnam.wing.semispan = Tecnam.wing.span/2;
Tecnam.wing.rootChord = 1.50;
Tecnam.wing.tipChord = Tecnam.wing.rootChord * Tecnam.wing.taperRatio;
Tecnam.wing.MAC = 2/3*Tecnam.wing.rootChord*(1 + Tecnam.wing.taperRatio + Tecnam.wing.taperRatio^2)/(1 + Tecnam.wing.taperRatio);
Tecnam.wing.surface = 2*(Tecnam.wing.semispan*Tecnam.wing.rootChord*(1 + Tecnam.wing.taperRatio)/2);
Tecnam.wing.surfaceProjected = Tecnam.wing.surface*cosd(Tecnam.wing.dihedralAngle);

Tecnam.wing.semispanwiseDiscr = 20; 
Tecnam.wing.chordwiseDiscr = 20;

Tecnam.wing.spanwiseSectors = 2 * Tecnam.wing.semispanwiseDiscr;
Tecnam.wing.delta_b = Tecnam.wing.span / Tecnam.wing.spanwiseSectors;
   
Tecnam.wing.LEpositionX = 0; 
Tecnam.wing.LEpositionY = 0; 
Tecnam.wing.LEpositionZ = 0;

Tecnam.wing.rotationangleX = 0;  % [deg]
% Tecnam.wing.rotationangleY = alpha;  defined as sequence for the Weissinger function
Tecnam.wing.rotationangleZ = 0;  % [deg]

% Tail structure data
Tecnam.tail.ID = 2;

Tecnam.tail.dihedralAngle = 0;  % [deg]
Tecnam.tail.sweepAngle = 0;  % [deg]

Tecnam.tail.aspectRatio = 3.49; 
Tecnam.tail.taperRatio = 1;
Tecnam.tail.span = 2.90; 
Tecnam.tail.semispan = Tecnam.tail.span/2;
Tecnam.tail.rootChord = 0.83;
Tecnam.tail.tipChord = Tecnam.tail.rootChord * Tecnam.tail.taperRatio;
Tecnam.tail.MAC = 2/3*Tecnam.tail.rootChord*(1 + Tecnam.tail.taperRatio + Tecnam.tail.taperRatio^2)/(1 + Tecnam.tail.taperRatio);
Tecnam.tail.surface = 2*(Tecnam.tail.semispan*Tecnam.tail.rootChord*(1 + Tecnam.tail.taperRatio)/2);
Tecnam.tail.surfaceProjected = Tecnam.tail.surface*cosd(Tecnam.tail.dihedralAngle);

Tecnam.tail.semispanwiseDiscr = 10;
Tecnam.tail.chordwiseDiscr = 10;

Tecnam.tail.spanwiseSectors = 2 * Tecnam.tail.semispanwiseDiscr;
Tecnam.tail.delta_b = Tecnam.tail.span / Tecnam.tail.spanwiseSectors;

Tecnam.tail.LEpositionX = 3.68;  % Positive to have the tail behind the wing
Tecnam.tail.LEpositionY = 0; 
Tecnam.tail.LEpositionZ = 0.28;  % Negative to have the tail below the wing

Tecnam.tail.rotationangleX = 0;  % [deg]
Tecnam.tail.rotationangleY = -2;  % [deg] (tail incidence)
Tecnam.tail.rotationangleZ = 0;  % [deg]

Tecnam.CL_3D = zeros(length(alpha_seq), 1);
Tecnam.CD_ind_3D = zeros(length(alpha_seq), 1);
for i = 1:length(alpha_seq)
    Tecnam.wing.rotationangleY = alpha_seq(i);
    [Gamma, CL, CD_ind, alpha_induced, elapsedTimeWeissinger] = fun_Weissinger(Tecnam.wing, Tecnam.tail, U_inf_modulus, U_inf, rho);
    
    % Data for polar plot (for AoA given range)
    Tecnam.CL_3D(i, 1) = CL.tot3D;
    Tecnam.CD_ind_3D(i, 1) = CD_ind.tot3D;
    
    % Data for circulation and CL distribution plot (at specific AoA)
    if Tecnam.wing.rotationangleY <= 3
       Tecnam.Wing_Grid_Gamma = Gamma.Wing_Grid_Gamma;
       Tecnam.Tail_Grid_Gamma = Gamma.Tail_Grid_Gamma;

       Tecnam.CL_2D_wing = CL.wing2D;
       Tecnam.CL_2D_tail = CL.tail2D;

       Tecnam.alpha_induced = alpha_induced;
    end
    
    % Elapsed time at last Weissinger function calls
    if Tecnam.wing.rotationangleY == alpha_seq(end)
       elapsedTimeWeissinger_last = elapsedTimeWeissinger;
    end
end


%% Plots section
% Circulation distribution plot
figure('Name','Circulation distribution',NumberTitle='off');
subplot(2, 2, 1); 
pcolor(Cessna.Tail_Grid_Gamma);
title('Cessna 172-SH tail circulation (\Gamma_{t})');
shading interp; 
colorbar;

subplot(2, 2, 3); 
pcolor(Cessna.Wing_Grid_Gamma); 
title('Cessna 172-SH wing circulation (\Gamma_{w})');
shading interp; 
colorbar;

subplot(2, 2, 2); 
pcolor(Tecnam.Tail_Grid_Gamma);
title('Tecnam P2002-JF tail circulation (\Gamma_{t})');
shading interp; 
colorbar;

subplot(2, 2, 4); 
pcolor(Tecnam.Wing_Grid_Gamma); 
title('Tecnam P2002-JF wing circulation (\Gamma_{w})');
shading interp; 
colorbar;


% Lift coefficient plot
Cessna.span_centers_wing = (-Cessna.wing.span/2 + Cessna.wing.delta_b/2):Cessna.wing.delta_b:(Cessna.wing.span/2 - Cessna.wing.delta_b/2);
Cessna.span_full_wing = linspace(-Cessna.wing.span/2, Cessna.wing.span/2, Cessna.wing.spanwiseSectors);

Cessna.span_centers_tail = (-Cessna.tail.span/2 + Cessna.tail.delta_b/2):Cessna.tail.delta_b:(Cessna.tail.span/2 - Cessna.tail.delta_b/2);
Cessna.span_full_tail = linspace(-Cessna.tail.span/2, Cessna.tail.span/2, Cessna.tail.spanwiseSectors);

Tecnam.span_centers_wing = (-Tecnam.wing.span/2 + Tecnam.wing.delta_b/2):Tecnam.wing.delta_b:(Tecnam.wing.span/2 - Tecnam.wing.delta_b/2);
Tecnam.span_full_wing = linspace(-Tecnam.wing.span/2, Tecnam.wing.span/2, Tecnam.wing.spanwiseSectors);

Tecnam.span_centers_tail = (-Tecnam.tail.span/2 + Tecnam.tail.delta_b/2):Tecnam.tail.delta_b:(Tecnam.tail.span/2 - Tecnam.tail.delta_b/2);
Tecnam.span_full_tail = linspace(-Tecnam.tail.span/2, Tecnam.tail.span/2, Tecnam.tail.spanwiseSectors);

figure('Name','Lift coefficient distribution',NumberTitle='off');
subplot(2, 2, 1);
plot(Cessna.span_centers_wing, Cessna.CL_2D_wing, '^g', LineWidth=2.0);
hold on
plot(Tecnam.span_centers_wing, Tecnam.CL_2D_wing, 'om', LineWidth=2.0);
title('Wing C_{l} at sector centers');
xlabel('wing span [m]'); 
ylabel('wing C_{l}');
legend('Cessna 172-SH', 'Tecnam P2002-JF',Location='south');
axis padded
grid on;

subplot(2, 2, 3);
plot(Cessna.span_full_wing, Cessna.CL_2D_wing, 'g', LineWidth=2.0);
hold on
plot(Tecnam.span_full_wing, Tecnam.CL_2D_wing, 'm', LineWidth=2.0);
title('Spanwise distribution interpolated');
xlabel('wing span [m]'); 
ylabel('wing C_{l}');
legend('Cessna 172-SH', 'Tecnam P2002-JF',Location='south');
axis padded
grid on;

subplot(2, 2, 2);
plot(Cessna.span_centers_tail, Cessna.CL_2D_tail, '^g', LineWidth=2.0);
hold on
plot(Tecnam.span_centers_tail, Tecnam.CL_2D_tail, 'om', LineWidth=2.0);
title('Tail C_{l} at sector centers');
xlabel('tail span [m]'); 
ylabel('tail C_{l}');
legend('Cessna 172-SH', 'Tecnam P2002-JF',Location='north');
axis padded
grid on;

subplot(2, 2, 4);
plot(Cessna.span_full_tail, Cessna.CL_2D_tail, 'g', LineWidth=2.0);
hold on
plot(Tecnam.span_full_tail, Tecnam.CL_2D_tail, 'm', LineWidth=2.0);
title('Spanwise distribution interpolated');
xlabel('tail span [m]'); 
ylabel('tail C_{l}');
legend('Cessna 172-SH', 'Tecnam P2002-JF',Location='north');
axis padded
grid on;


% Induced angle plot
figure('Name','Induced angle distribution',NumberTitle='off');
subplot(2, 2, 1);
plot(Cessna.span_centers_wing, rad2deg(Cessna.alpha_induced(1, 1:Cessna.wing.spanwiseSectors)), '^g', LineWidth=2.0);
hold on
plot(Tecnam.span_centers_wing, rad2deg(Tecnam.alpha_induced(1, 1:Tecnam.wing.spanwiseSectors)), 'om', LineWidth=2.0);
title('Induced angle on wing');
xlabel('wing span [m]'); 
ylabel('wing \alpha_{ind} [deg]');
legend('Cessna 172-SH', 'Tecnam P2002-JF',Location='south');
axis padded
grid on;

subplot(2, 2, 3);
plot(Cessna.span_full_wing, rad2deg(Cessna.alpha_induced(1, 1:Cessna.wing.spanwiseSectors)), 'g', LineWidth=2.0);
hold on
plot(Tecnam.span_full_wing, rad2deg(Tecnam.alpha_induced(1, 1:Tecnam.wing.spanwiseSectors)), 'm', LineWidth=2.0);
title('Induced angle distribution interpolated');
xlabel('wing span [m]'); 
ylabel('wing \alpha_{ind} [deg]');
legend('Cessna 172-SH', 'Tecnam P2002-JF',Location='south');
axis padded
grid on;

subplot(2, 2, 2);
plot(Cessna.span_centers_tail, rad2deg(Cessna.alpha_induced(1, (Cessna.wing.spanwiseSectors + 1) : (Cessna.wing.spanwiseSectors + Cessna.tail.spanwiseSectors))), '^g', LineWidth=2.0);
hold on
plot(Tecnam.span_centers_tail, rad2deg(Tecnam.alpha_induced(1, (Tecnam.wing.spanwiseSectors + 1) : (Tecnam.wing.spanwiseSectors + Tecnam.tail.spanwiseSectors))), 'om', LineWidth=2.0);
title('Induced angle on tail');
xlabel('tail span [m]'); 
ylabel('tail \alpha_{ind} [deg]');
legend('Cessna 172-SH', 'Tecnam P2002-JF',Location='north');
axis padded
grid on;

subplot(2, 2, 4);
plot(Cessna.span_full_tail, rad2deg(Cessna.alpha_induced(1, (Cessna.wing.spanwiseSectors + 1) : (Cessna.wing.spanwiseSectors + Cessna.tail.spanwiseSectors))), 'g', LineWidth=2.0);
hold on
plot(Tecnam.span_full_tail, rad2deg(Tecnam.alpha_induced(1, (Tecnam.wing.spanwiseSectors + 1) : (Tecnam.wing.spanwiseSectors + Tecnam.tail.spanwiseSectors))), 'm', LineWidth=2.0);
title('Induced angle distribution interpolated');
xlabel('tail span [m]'); 
ylabel('tail \alpha_{ind} [deg]');
legend('Cessna 172-SH', 'Tecnam P2002-JF',Location='north');
axis padded
grid on;


% Polar plot
figure('Name','Polar plot',NumberTitle='off');
plot(Cessna.CD_ind_3D, Cessna.CL_3D, '-^g', LineWidth=2.0);
hold on
plot(Tecnam.CD_ind_3D, Tecnam.CL_3D, '-om', LineWidth=2.0);
title('Lift-Drag polar curve');
xlabel('C_{D, induced} (induced Drag coefficient)');
ylabel('C_{L} (Lift coefficient)');
legend(['Cessna 172 SkyHawk (tail incidence = ' num2str(Cessna.tail.rotationangleY, '%.1f') ' [deg])'],['Tecnam P2002-JF (tail incidence = ' num2str(Tecnam.tail.rotationangleY, '%.1f') ' [deg])'], Location='southeast', FontSize=15);
axis padded
grid on


% Elapsed time for computation
fprintf('Elapsed time first call Weissinger function: %.4f [s]\n', elapsedTimeWeissinger_first);
fprintf('Elapsed time last call Weissinger function: %.4f [s]\n\n', elapsedTimeWeissinger_last);

Cessna.E = Cessna.CL_3D./Cessna.CD_ind_3D;
Tecnam.E = Tecnam.CL_3D./Tecnam.CD_ind_3D;