%CU Spring 2024 Aerospace Dynamics Project
clear;close all; clc

%==========================================================================
%                        ** OPERATIONAL PARMS. ** 
%Blade Params.
m_b = 6.4006;     %Mass of Blade (kg)
k_b = 293.0275e6; %Spring Constant of Blade (Nm)
n   = 16;         %Number of Blades
i   = 1:n;        %Individual Blade Designation

%Engine Parms.
omega = 262.85;            %Engine Rotational Speed (rad/s)
RPM   = omega*60 / (2*pi); %RPM of engine fan, (rev/min)

%External Force Parms.
g   = 9.81;                       %Acceleration due to gravity (m/s^2)
D_b = Blade_Drag(RPM, 0.54879); %Drag force on an individual blade, N

%Blade Defect Params
d_b = 12; %Defect Blade Number
blade_defect_percent = 1.0;


%==========================================================================
%                       ** MATRIX ASSEMBLY ** 
M_blade              = (m_b)*eye(2*n);
M_blade(d_b,d_b)     = M_blade(d_b,d_b)*blade_defect_percent;
M_blade(n+d_b,n+d_b) = M_blade(n+d_b,n+d_b)*blade_defect_percent;

K_blade = k_b*[eye(n) zeros(n);zeros(n) zeros(n)];

F_blade = D_b*[zeros(n,1);ones(n,1)];
F_blade(n+d_b,1) = F_blade(n+d_b,1)*(blade_defect_percent^4);

%Mass Matrix
M = [0.5*trace(M_blade)*eye(2) [cos((2*pi*(i-1))/n) -sin((2*pi*(i-1))/n);sin((2*pi*(i-1))/n) cos((2*pi*(i-1))/n)]*M_blade;M_blade*[cos((2*pi*(i-1))/n)' sin((2*pi*(i-1))/n)';-sin((2*pi*(i-1))/n)' cos((2*pi*(i-1))/n)'] M_blade];

%Stiffness Matrix
K_spring = [trace(K_blade)*eye(2) [cos((2*pi*(i-1))/n) -sin((2*pi*(i-1))/n);sin((2*pi*(i-1))/n) cos((2*pi*(i-1))/n)]*K_blade;K_blade*[cos((2*pi*(i-1))/n)' sin((2*pi*(i-1))/n)';-sin((2*pi*(i-1))/n)' cos((2*pi*(i-1))/n)'] K_blade];
K_centripetal = (omega^2)*M;
K = K_spring+K_centripetal;

%Force Matrix
F = [eye(2) [cos((2*pi*(i-1))/n) -sin((2*pi*(i-1))/n);sin((2*pi*(i-1))/n) cos((2*pi*(i-1))/n)];[cos((2*pi*(i-1))/n)' sin((2*pi*(i-1))/n)';-sin((2*pi*(i-1))/n)' cos((2*pi*(i-1))/n)'] eye(2*n)]*[zeros(2,1);F_blade];


%==========================================================================
%                         ** ANALYSIS ** 
[EVec, Eval, NatFreq, mu, gamma] = MDOF_Analysis(M,K);

NatFreq








