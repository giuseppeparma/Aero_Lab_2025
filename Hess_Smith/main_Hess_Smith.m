clear
close all
clc

%% Hess-Smith method initialization
% addpath mat_functions

U_inf = 1;                  % Far-field velocity [m/s]
alpha = 1;                  % Angle of attack [deg]
U_inf_x = U_inf * cos(deg2rad(alpha));
U_inf_y = U_inf * sin(deg2rad(alpha));
U_inf = [U_inf_x; U_inf_y];

Chord = 1;
NPanels = 160;

naca='0008';
Re = 1e6;
iter=250;

LE_X_Position = 0;
LE_Y_Position = 0;

% Profile creation from XFoil
[x,y] = createXfoil('profile',naca, NPanels, Chord);

geo.x = x;
geo.y = y;

figure('Name','Airfoil panelization',NumberTitle='off');
plot(x, y, '-o');
title('Profile');
axis equal
grid on

% Create discretisation & initialisation of problem matrices
[centers,normals,tangents,extrema_1,extrema_2,alpha_panel,lengths,L2G_TransfMatrix,G2L_TransfMatrix] = createPanels(geo);
        
NCols = NPanels + 1;
NRows = NCols;
A = zeros(NRows, NCols);     % system coefficients
B = zeros(NRows, 1);         % known terms

% A matrix
% non-penetration equations
for i = 1:NPanels
    local_center = centers(i, :)';
    local_normal = normals(i, :)';

    for j = 1:NPanels
        local_extreme_1 = extrema_1(j, :)';
        local_extreme_2 = extrema_2(j, :)';

        local_L2G_TransfMatrix = squeeze(L2G_TransfMatrix(j, :, :));
        local_G2L_TransfMatrix = squeeze(G2L_TransfMatrix(j, :, :));

        A(i, j) = dot(uSource(local_center, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), local_normal);

        A(i, sum(NPanels)+1) = A(i, sum(NPanels)+1) + dot(uVortex(local_center, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), local_normal);
    end
end

% Kutta condition equation
first_centers = centers(1, :)';
first_tangents = tangents(1, :)';

last_centers = centers(end, :)';
last_tangents = tangents(end, :)';

last_a = 0;
for j = 1:NPanels
    local_extreme_1 = extrema_1(j, :)';
    local_extreme_2 = extrema_2(j, :)';
    local_L2G_TransfMatrix = squeeze(L2G_TransfMatrix(j, :, :));
    local_G2L_TransfMatrix = squeeze(G2L_TransfMatrix(j, :, :));

    a = dot(uSource(first_centers, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), first_tangents);
    last_a = last_a + dot(uVortex(first_centers, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), first_tangents);

    a = a + dot(uSource(last_centers, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), last_tangents);
    last_a = last_a + dot(uVortex(last_centers, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix), last_tangents);

    A(NPanels + 1, j) = a;
end

A(NPanels + 1, NPanels + 1) = last_a;

% B matrix
for j = 1:NPanels
    local_normal = normals(j, :)';
    B(j) = - dot(U_inf, local_normal);
end

B(NPanels + 1) = - dot(U_inf, (first_tangents + last_tangents));

% Solution of the linear system
solution = linsolve(A, B);


% Velocity field computation
% Mesh
x_p = linspace(-0.2, 1.2, 20);
y_p = linspace(-0.3, 0.3, 20);
diff = zeros(NPanels);

U_tot = zeros(20,20);
V_tot = zeros(20,20);

figure('Name','Velocity field',NumberTitle='off');
plot(x, y, '-o');
hold on
title('Velocity field around the profile');
axis equal
grid on

for p = 1:20
    for q = 1:20
        point = [x_p(p); y_p(q)]; % point of velocity computation

        % Check if the point is inside or on the contour of the profile
        [in,on] = inpolygon(x_p(p), y_p(q), x, y);
        if in || on
            continue
        end

        source_induced_tot = [0; 0];
        vortex_induced_tot = [0; 0];
        source_panel = [0; 0];
        vortex_panel = [0; 0];
            for k = 1:NPanels
                Temp_L2G_TransfMatrix = squeeze(L2G_TransfMatrix(k, :, :));
                Temp_G2L_TransfMatrix = squeeze(G2L_TransfMatrix(k, :, :));

                source_panel = solution(k)*uSource(point,extrema_1(k,:)',extrema_2(k,:)',Temp_L2G_TransfMatrix,Temp_G2L_TransfMatrix);
                source_induced_tot = source_induced_tot + source_panel;

                vortex_panel = solution(NPanels + 1)*uVortex(point,extrema_1(k,:)',extrema_2(k,:)',Temp_L2G_TransfMatrix,Temp_G2L_TransfMatrix);
                vortex_induced_tot = vortex_induced_tot + vortex_panel;
            end

        U_tot(p,q) = U_inf_x + source_induced_tot(1) + vortex_induced_tot(1);
        V_tot(p,q) = U_inf_y + source_induced_tot(2) + vortex_induced_tot(2);
    end
end

for i = 1:20
    for j = 1:20
        scale = 0.09;
        quiver(x_p(i), y_p(j), U_tot(i,j)*scale, V_tot(i,j)*scale);
    end
end

%% Hess Smith method validation

% Cp computation
% Hess-Smith
Cp_HS = pressureCoeff(NPanels, centers, extrema_1, extrema_2, tangents, L2G_TransfMatrix, G2L_TransfMatrix, solution, U_inf);
% XFoil 
[x_Cp, Cp_XF] = createXfoil('CP',naca,NPanels,Chord,alpha);

figure('Name','Pressure coefficient',NumberTitle='off');
plot(centers(:,1), Cp_HS, 'm');
hold on
plot(x_Cp, Cp_XF, 'g');
set(gca, 'YDir', 'reverse');
title ('Pressure coefficient comparison between Hess-Smith implementation and XFoil results');
xlabel ('Chord position');
ylabel ('Cp');
legend('Hess-Smith implemented Cp', 'XFoil computed Cp', Location='best'); 
grid on

% Error between Hess-Smith and XFoil pressure coefficient
Cp_err = errCp(NPanels,centers,x_Cp,Cp_HS,Cp_XF);
figure('Name', 'Error evaluation',NumberTitle='off');
plot(x_Cp,Cp_err,'-*r');
title('Cp_{HS} - Cp_{XF}');
xlabel('Chord position by panels extrema');
ylabel('Error');
axis padded
grid on


% Cl computation
% Hess Smith
Cl = 0;
normal2U_inf = [sin(deg2rad(alpha));-cos(deg2rad(alpha))]; % normal direction to the velocity
for i = 1:NPanels
    Cl = Cl + Cp_HS(i)*(lengths(i)/Chord)*dot(normals(i,:),normal2U_inf);
end

fprintf('\nComputed lift coefficient CL = %f\n', Cl);

% Xfoil
[coeffs,]=createXfoil('Coeffs',naca,NPanels,Chord,alpha);
xfoilCL=coeffs.cl;
fprintf('\nXfoil lift coefficient CL = %f\n\n',xfoilCL);



% Cm computation
% Hess Smith
Cm = 0;
z = [0;0;1];

for i = 1:NPanels
    r = [centers(i,1) - Chord/4; centers(i,2); 0];  % centroyd - chord/4 distance
    normal3d = [normals(i,1); normals(i,2); 0];

    Cm = Cm - Cp_HS(i)*(lengths(i)/Chord)*dot(cross(r,normal3d),z); 
end

fprintf('\nComputed moment coefficient Cm = %f\n', Cm);

% Xfoil
% Beware of the conventions: the sign is different due to different choice of axes
xfoilCM=coeffs.cm; 
fprintf('\nXfoil moment coeffficient = %f\n',xfoilCM);




%% Comparison between NACA0008 and NACA1410
% Creation  of the airfoils: the first one has already been defined, now
% define the comparison one

naca1='1410';

LE_X_Position1 = 0;
LE_Y_Position1 = 0;

% Profile creation from XFoil
[x1,y1] = createXfoil('profile',naca1, NPanels, Chord);

geo1.x = x1;
geo1.y = y1;

% Definition of sequence Reynolds numbers
ReSeq=[1,5,10]*1e6; 


% Transition laminar-to-turbulent

% Initialisation of the variables
x_tr=zeros(length(ReSeq));
x1_tr=zeros(length(ReSeq));

% Computation of the values for each Reynolds for both airfoils
for k=1:length(ReSeq)
    coeffs_tr=createXfoil('Coeffsvisc',naca,NPanels,Chord,alpha,ReSeq(k),iter);
    x_tr(k)=coeffs_tr.top_xtr;
    
    coeffs1_tr=createXfoil('Coeffsvisc',naca1,NPanels,Chord,alpha,ReSeq(k),iter);
    x1_tr(k)=coeffs1_tr.top_xtr;
end

% The displayed points are chord positions, therefore plot them on the x-axis
y_tr=zeros(length(x_tr),1);

% Comparison plots with different turbulence levels
figure('Name','Airfoil comparison: chord position of transition point',NumberTitle='off')
subplot(1,2,1)
plot(x, y,'b','LineWidth',2);
hold on
for k=1:length(ReSeq)
plot(x_tr(k,1),y_tr,'*','LineWidth',1.5);
end
legend('airfoil',['Re = ' num2str(ReSeq(1))],['Re = ' num2str(ReSeq(2))],['Re = ' num2str(ReSeq(3))],'Location','best')
title('NACA0008')
daspect([1 1 1])
grid on

subplot(1,2,2)
plot(x1, y1,'b','LineWidth',2);
hold on
for k=1:length(ReSeq)
plot(x1_tr(k,1),y_tr,'*','LineWidth',1.5);
end
legend('airfoil',['Re = ' num2str(ReSeq(1))],['Re = ' num2str(ReSeq(2))],['Re = ' num2str(ReSeq(3))],'Location','best')
title('NACA1410')
daspect([1 1 1])
grid on



% Separation

figure('Name','Airfoil comparison: separation',NumberTitle='off')
% Computation of friction coefficient for the required series of Reynolds numbers
% Each subplot showcases both suction and pressure sides
for n=1:length(ReSeq)
    [CfTop,CfBot]=createXfoil('CF',naca,NPanels,Chord,alpha,ReSeq(n),iter);
    subplot(2,3,n)
    plot(CfTop(:,1),CfTop(:,2),'g');
    hold on
    plot(CfBot(:,1),CfBot(:,2),'m');
    hold on
    title(['NACA0008, Re = ' num2str(ReSeq(n))]);
    axis padded
    xlabel x
    ylabel Cf
    grid on
    legend('suction side','pressure side','Location','ne')
end
hold on

% The computations are carried out for both studied profiles, and then plotted together
for n=1:length(ReSeq)
    [Cf1Top,Cf1Bot]=createXfoil('CF',naca1,NPanels,Chord,alpha,ReSeq(n),iter);
    subplot(2,3,n+3)
    plot(Cf1Top(:,1),Cf1Top(:,2),'g');
    hold on
    plot(Cf1Bot(:,1),Cf1Bot(:,2),'m');
    hold on
    title(['NACA1410, Re = ' num2str(ReSeq(n))]);
    axis padded
    xlabel x
    ylabel Cf
    grid on 
    legend('suction side','pressure side','Location','ne')
end


