function [D] = Blade_Drag(RPM, air_density)
%{
This function calculates the total drag on a single engine blade. 
INPUTS
    RPM = Engine rotational rate in rev/minute
    density = air density in kg/m^3
OUTPUTS
    D = Blade Drag in N
%}

    %Blade Geometry and Properties:
    c  = 4*(1/12)*(1/3.28); %Length averaged chord, [m]
    l  = 0.925*(3.4044/2);  %Blade length corrected for hub radius, [m]
    CD = 1.0;               %Assumed Drag Coefficient
    
    %Engine Parms.
    omega = RPM*2*pi/60;
    
    %Calculations
    S = l*c;     %Blade planform area, [m^2]
    D = (1/6) * l^3 * CD * air_density * omega^2 * S; %Single Blade Drag
end