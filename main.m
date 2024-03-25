%CU Spring 2024 Aerospace Dynamics Project
clear;close all; clc

%Rotor Params
m_r = 1; %Mass of Rotor (kg)
k_r = 1; %Spring Constant of Rotor

%Blade Params
m_b = 1; %Mass of Blade (kg)
k_b = 1; %Spring Constant of Blade

n = 12; %Number of Blades
i = 1:n; %Individual Blade Designation

%External Forces Params
omega = 200; %Rotational Speed (rad/s)
g = 9.81; %Acceleration due to gravity (m/s^2)
D_b = 15; %Drag Force (N)

%Alternative Eye Matrix
alt_eye = zeros(2*n);
for e = 1:size(alt_eye,1) 
    if mod(e,2) == 1 
        alt_eye(e,e) = 1;
    end
end

%Alternative One Matrix
alt_ones = zeros(2*n,1);
for e = 1:size(alt_ones,1) 
    if mod(e,2) == 0 
        alt_ones(e,1) = 1;
    end
end

%Mass Matrix
M = [(m_r+n*m_b)*eye(2) m_b*[cos((2*pi*(i-1))/n) -sin((2*pi*(i-1))/n);sin((2*pi*(i-1))/n) cos((2*pi*(i-1))/n)];m_b*[cos((2*pi*(i-1))/n)' sin((2*pi*(i-1))/n)';-sin((2*pi*(i-1))/n)' cos((2*pi*(i-1))/n)'] (m_b)*eye(2*n)];

%Spring Matrix
K_spring = [(k_r)*eye(2) zeros(2,2*n);zeros(2*n,2) k_b*alt_eye];
K_centripetal = -M*omega^2;
K = K_spring+K_centripetal;

%External Force
F_g = [eye(2) [cos((2*pi*(i-1))/n) -sin((2*pi*(i-1))/n);sin((2*pi*(i-1))/n) cos((2*pi*(i-1))/n)];[cos((2*pi*(i-1))/n)' sin((2*pi*(i-1))/n)';-sin((2*pi*(i-1))/n)' cos((2*pi*(i-1))/n)'] eye(2*n)]*g*[m_r*[0 -1]';m_b*[-sin((2*pi*(i-1))/n) -cos((2*pi*(i-1))/n)]'];
F_aero = [eye(2) [cos((2*pi*(i-1))/n) -sin((2*pi*(i-1))/n);sin((2*pi*(i-1))/n) cos((2*pi*(i-1))/n)];[cos((2*pi*(i-1))/n)' sin((2*pi*(i-1))/n)';-sin((2*pi*(i-1))/n)' cos((2*pi*(i-1))/n)'] eye(2*n)]*[[0 0]';D_b*alt_ones];
