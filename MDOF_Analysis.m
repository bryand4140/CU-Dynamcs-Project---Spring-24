function [EVec, Eval, NatFreq, mu, gamma] = MDOF_Analysis(M,K)
%This function computes the eigenvectors and eigenvalues for a MDOF
%system which is described by the mass matrix M and the stiffness matrix K.
%{
INPUTS:
    M = The nxn mass matrix
    K = The nxn stiffness matrix

OUTPUTS:
EVec    = A nxn matrix whose columns are the eigenvectors for the system.
Eval    = A nx1 column vector whose elements are the eigenvalues for 
          the system.
NatFreq = A nx1 column vector whose elements are the natural frequencies 
          of the system. 
mu      = A nx1 column vector whose elements are the normalized modal
          masses.
gamma   = A nx1 column vector whose elements are normalized the model 
          stiffnesses.
%}
    
    
    %Calculate the Eigenvectors:
    [EVec, D] = eig(K, M);
    
    %Calculate the eigenvalues and natural frequencies
    Eval    = diag(D);
    NatFreq = sqrt(Eval); %Rad/s
    
    %Calculate the geeralized mass and stiffness values:
    n = size(EVec,1); %Determine the size of the system (number of DOF)
    
    %Preallocate:
    mu    = zeros(n,1);
    gamma = zeros(n,1);
    
    for i = 1:n
        mu(i)    = transpose(EVec(:,i))*M*EVec(:,i);
        gamma(i) = transpose(EVec(:,i))*K*EVec(:,i);
    end

end %End of MDOF_Analysis 