function [D] = Blade_Aero_Force(omega,air_density,l,CA)
%{
This function calculates the total drag on a single engine blade. 
INPUTS
    omega = Engine rotational rate in rad/s
    density = air density in kg/m^3
OUTPUTS
    D = Blade Drag in N
%}

    %Blade Geometry and Properties:
    c  = 4*(1/12)*(1/3.28); %Length averaged chord, [m]
    
    %Calculations
    D = (1/6) * air_density * omega^2 * c *CA * l^4;
end