clear;close all;clc
format shorteng

%Code Description: This code calculates the frequency response of a turbofan
%engine. The force on the engine is calculated by considering centripetal
%force and aerodynamic force on the blades. This translates to an external
%force in the x and y direction which are the two degrees of freedom. 

%Define system paramters
mb  = 6.4005;  %Mass of a single blade, [kg]
me  = 5.3e3;   %Mass of the engine, [kg]
kx  = 500;     %Spring constant for the engine (x-dir), [N/m]
ky  = 500;     %Spring constant for the engine (y-dir), [N/m]
w   = 250.0;   %Rotational speed of the rotor, [rad/s]
rho = 0.01841; %Air density, [kg/m^3]
nb  = 8;       %Number of blades
t   = 0.02;    %Thickness of the blade, [m]
c   = 0.1016;  %Chord length of the blade, [m]
L   = 1.5745;  %Length of the blade, [m]
CD  = 0.1;     %Drag coefficient

%Mass and stiffness matricies:
M = [me, 0 ; 0, me];
K = [kx 0 ; 0 ky];

%Initial Conditions:
x0 = [0;0]; %Initial x and y positions
v0 = [0;0]; %Initial xdot and ydot velocities

%Perform Eign Analysis
[EVec, Eval, NatFreq, mu, gamma] = MDOF_Analysis(M,K);

%Create arrays for the position angles, mass of the blades,and length:
theta_blade = (2*pi/nb)*(180/pi)*(0:1:nb);
m_blades    = mb * ones(1,nb);
L_blades    = L  * ones(1,nb);

%Modify Blade 4:
%Reduce mass and length by 25% to simulate failure
m_blades(4) = mb*0.75;
L_blades(4) = L*0.75;

F_ext_blades = zeros(2,nb);
for i = 1:nb
    %Form transformation matrix for blade i
    R_BN = [cosd(theta_blade(i)), -sind(theta_blade(i));sind(theta_blade(i)), cosd(theta_blade(i))];

    %Calculate aero and centripital force for blade i
    F_aero = Blade_Aero_Force(w,rho,L_blades(i),CD);
    F_cent = Blade_Cent_Force(L_blades(i),m_blades(i),w);

    %Calculate the total force vector in blade i
    F_ext_i = R_BN*[F_cent; -1*F_aero];

    %Assign to external force storage matrix
    F_ext_blades(1,i) = F_ext_i(1);
    F_ext_blades(2,i) = F_ext_i(2);
end

%Extract elements and sum in x and y direction:
Fx = sum(F_ext_blades(1,:)); %Total force in x-direction, [N]
Fy = sum(F_ext_blades(2,:)); %Total force in y-direction, [N]

%Define the total time span
tspan = 0:0.1:20;

%Calculate Particular Solutions:
xp = zeros(1,length(tspan));
yp = zeros(1,length(tspan));
for i = 1:length(tspan)
    xp(i) = (Fx*cosh(sqrt(-kx/me)*tspan(i)))/(me*w^2-kx) - (Fx*cos(tspan(i)*w))/(me*w^2-kx);
    yp(i) = (Fy*cosh(sqrt(-ky/me)*tspan(i)))/(me*w^2-kx) - (Fy*cos(tspan(i)*w))/(me*w^2-ky);
end

figure('Color','white')
subplot(2,1,1)
plot(tspan,xp*1e3,'b','LineWidth',1.5)
xlabel('Time, s')
ylabel('Response Amplitude., mm')
title('x-Direction')

subplot(2,1,2)
plot(tspan,xp*1e3,'r','LineWidth',1.5)
xlabel('Time, s')
ylabel('Response Amplitude, mm')
title('y-Direction')

function [Response] = calculateResponse(M, K, tspan, x0, v0, Fx, Fy, w)
    % This function calculates the dynamic response of a turbofan engine
    % to external forces in the x and y directions as a function of time and
    % returns a single matrix combining the time vector and displacement responses.
    %
    % INPUTS:
    % M      - The mass matrix
    % K      - The stiffness matrix
    % tspan  - The time span for the simulation
    % x0     - Initial displacement vector
    % v0     - Initial velocity vector
    % Fx     - Constant external force in the x direction
    % Fy     - Constant external force in the y direction
    % w      - Angular velocity
    %
    % OUTPUT:
    % Response - A matrix where the first column is the time vector, 
    %            the second column is the x displacement response over time, 
    %            and the third column is the y displacement response over time.

    % Initialization
    dt = tspan(2) - tspan(1); % Time step
    n = length(tspan); % Number of time steps
    X_t = zeros(2, n); % Displacement response over time
    X_t(:, 1) = x0; % Initial condition
    
    % Perform Eigen Analysis (to obtain system frequencies, damping can be added if necessary)
    [EVec, ~, NatFreq, ~, ~] = MDOF_Analysis(M, K);

    % Main simulation loop
    for i = 2:n
        t = tspan(i);
        integral_x = zeros(2, 1);
        integral_y = zeros(2, 1);
        for tau = tspan(1):dt:t-dt
            % Force contributions at each time step
            Fx_tau = Fx * cos(w * tau);
            Fy_tau = Fy * cos(w * tau);

            % Duhamel's integral approximation for each mode
            integral_x = integral_x + Fx_tau * sin(NatFreq * (t - tau)) * dt;
            integral_y = integral_y + Fy_tau * sin(NatFreq * (t - tau)) * dt;
        end
        
        % Calculate total displacement for this time step
        for mode = 1:2 % For each mode
            X_t(:, i) = X_t(:, i) + EVec(:, mode) .* (integral_x(mode) + integral_y(mode)) / (M(mode, mode) * NatFreq(mode)^2);
        end
    end
    
    % Concatenate the time vector with displacement responses to form the output matrix
    Response = [tspan.', X_t(1, :).', X_t(2, :).'];
end








