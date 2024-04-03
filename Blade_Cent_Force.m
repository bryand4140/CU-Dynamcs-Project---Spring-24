function [Fc] = Blade_Cent_Force(L,m,omega)
%{
This function computes the centripital force on a given blade.

INPUTS:
    L = length of the blade, [m]
    m = mass of the blade, [m]
    omega = angular velocity, [rad/s]
 OUTPUTS:
    Fc = Centripital Force, [N]
%}

Fc = (1/2) * L * m * omega^2;
end
