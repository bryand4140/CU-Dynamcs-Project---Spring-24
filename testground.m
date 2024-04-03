clear;close all;clc
format shorteng

%Code Description: This code calculates the frequency response of a turbofan
%engine. The force on the engine is calculated by considering centripetal
%force and aerodynamic force on the blades. This translates to an external
%force in the x and y direction which are the two degrees of freedom. 

%Define system paramters
mb  = 6.4005;  %Mass of a single blade, [kg]
me  = 5.3e3;   %Mass of the engine, [kg]
kex = 1e3;     %Spring constant for the engine (x-dir), [N/m]
key = 1e3;     %Spring constant for the engine (y-dir), [N/m]
w   = 250.0;   %Rotational speed of the rotor, [rad/s]
rho = 0.01841; %Air density, [kg/m^3]
nb  = 8;       %Number of blades
t   = 0.02;    %Thickness of the blade, [m]
c   = 0.1016;  %Chord length of the blade, [m]
L   = 1.5745;  %Length of the blade, [m]
CD  = 0.1;     %Drag coefficient

%Mass and stiffness matricies:
M = [me, 0 ; 0, me];
K = [kex 0 ; 0 key];

%Perform Eign Analysis
[EVec, Eval, NatFreq, mu, gamma] = MDOF_Analysis(M,K);


%Create arrays for the position angles, mass of the blades,and length:
theta_blade = (2*pi/nb)*(180/pi)*(0:1:nb);
m_blades    = mb * ones(1,nb);
L_blades    = L  * ones(1,nb);

%Modify Blade 4:
%Reduce mass and length by 25%
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
Fx = sum(F_ext_blades(1,:));
Fy = sum(F_ext_blades(2,:));

%Define the total time span
tspan = 0:0.1:5;

%Pre-allocate arrays for the forces in the x and y as a function of time:
Fx_t = zeros(1,length(tspan));
Fy_t = zeros(1,length(tspan));

for i = 1:length(tspan)
    Fx_t(i) = Fx*cos(w*tspan(i));
    Fy_t(i) = Fy*cos(w*tspan(i));
end


[x_response, y_response] = Response_Duhmel(M, K, tspan, Fx_t, Fy_t);



function [x_response, y_response] = Response_Duhmel(M, K, tspan, Fx_t, Fy_t)
    % M: Mass matrix
    % K: Stiffness matrix
    % tspan: Time span for the simulation
    % Fx_t: External force in the x direction as a function handle
    % Fy_t: External force in the y direction as a function handle
    
    % Compute modal analysis
    [EVec, Eval, ~, mu, ~] = MDOF_Analysis(M, K);
    NatFreq = sqrt(Eval); % Natural frequencies
    
    % Number of DOFs
    n = length(NatFreq);
    
    % Time vector
    t = linspace(tspan(1), tspan(2), length(tspan));
    
    % Initialize responses
    x_response = zeros(1, length(t));
    y_response = zeros(1, length(t));
    
    % Loop through each mode
    for i = 1:n
        % Modal mass for the i-th mode
        m_i = mu(i);
        
        % Natural frequency for the i-th mode
        omega_i = NatFreq(i);
        
        % Impulse response function for the i-th mode (no damping)
        h_i = @(tau) (1 / (m_i * omega_i)) * sin(omega_i * tau);
        
        % Duhamel Integral for each time step and each direction
        for k = 1:length(t)
            % Time up to current step
            tk = t(k);
            
            % Integration for x and y directions
            x_integral = integral(@(tau) h_i(tk - tau) .* Fx_t(tau), 0, tk, 'ArrayValued', true);
            y_integral = integral(@(tau) h_i(tk - tau) .* Fy_t(tau), 0, tk, 'ArrayValued', true);
            
            % Aggregate modal responses (scaled by eigenvectors)
            x_response(k) = x_response(k) + EVec(1,i) * x_integral;
            y_response(k) = y_response(k) + EVec(2,i) * y_integral;
        end
    end
end







