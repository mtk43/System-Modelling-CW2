function [Melem] = LocalMassMatrix(eID, msh, gq)
%Function which generates the local element mass matrix (Melem), of the linear
%reaction term
%%Takes in the element ID (eID), the 1D mesh (msh),and the Gaussian Quadrature scheme (GQ)
    %% Initialise necessary variables 
    N = 2;              % Define the GQ level
    Melem = zeros(2,2);   % Initialise an empty 2x2 LEM matrix
    
    %% Determine values
    % Loop through Gauss points
    for it = 1:N
        % Evaluate the basis function
        psi_0 = EvalBasis(0, gq.xipts(it));
        psi_1 = EvalBasis(1, gq.xipts(it));

        % Generate matrix with psi values in correct indices
        psi = [psi_0*psi_0, psi_0*psi_1; psi_1*psi_0, psi_1*psi_1];
        
        % Multiply by Gauss weight
        psi = psi * gq.gsw(it);
        
        % Adding onto the matrix
        Melem = Melem + psi;
        
    end
    Melem = Melem * msh.elem(eID).J;
end