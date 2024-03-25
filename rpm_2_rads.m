function [out] = rpm_2_rads(rpm)
%This function converts an input in RPM to an output in radians/second

out = rpm*(2*pi/60);

end