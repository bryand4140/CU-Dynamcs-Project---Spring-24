clear;close all; clc; format shorteng

t = 0.02;              %Chord averaged thickness, [m]
c = 4*(1/12)*(1/3.28); %Length averaged chord, [m]
l = 0.925*(3.4044/2);  %Blade length corrected for hub radius, [m]

%Material properties
rho_b = 2.0e3;   %Approx. carbon composite density, [kg/m^3]
E     = 227.0e9; %Approx. Young's modulus for carbon fiber in tension, [N/m^2]

%Calculations
A = t*c;     %Length averaged cross sectional area, [m^2]
V = t*c*l;   %Approximate blade volume, [m^3]
S = l*c;     %Blade planform area, [m^2]
m = rho_b*V;   %Mass of a single blade, [kg]
k_t = E*A/l; %Spring constant for a single blade, [N/m]


%Print Statements
fprintf('Mass of a single blade = %0.4f kg\n',m)
fprintf('Blade tensile spring constant = %0.4f MN/m\n',k_t/1e6)