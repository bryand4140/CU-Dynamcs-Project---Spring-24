function [out] = Norm_NatFreq_MaxVal(x)
%This function normalizez all values (frequencies) in an array x by the
%maxximum value in the array.

    max_x = max(x);
    out   = zeros(1,length(x));
    
    for i = 1:length(x)
        out(i) = x(i)/max_x;
    end
end
