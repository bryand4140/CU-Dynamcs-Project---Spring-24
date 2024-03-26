function [Fc] = Centripital_Blade_Force(RPM)
%{
This function calculates the total centripital force on a single blade for 
a given RPM.

INPUTS
    RPM = Engine rotational rate in revolutions/minute
OUTPUTS
    Fc = Centripital Force in N
%}

    %Blade Geometry and Properties:
    L = 0.925*(3.4044/2);  %Blade length corrected for hub radius, [m]
    m = 6.4005; %kg

    %Engine Parms.
    omega = rpm_2_rads(RPM);

    Fc = 0.5 * L* m * omega^2;
end
