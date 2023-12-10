function [Kelem, diffusion_term] = LocalStiffnessMatrix(eID, msh, D, lambda, gq)
%%This function uses constant coefficient values (should definitely be
%%generalised later)
% Extract J value
J = msh.elem(eID).J;
%% Solving diffusion term
diffusion_term = zeros(2,2);
N = 2;
% Loop over the Gauss xipts and weights
for i = 1:N
    
    % Evaluate the basis functions at the current Gauss point
    psi_deriv_0 = EvalBasisGrad(0);
    psi_deriv_1 = EvalBasisGrad(1);
    
    % Form the local element matrix for the diffusion term
    psi_deriv = [psi_deriv_0*psi_deriv_0, psi_deriv_0*psi_deriv_1; psi_deriv_1*psi_deriv_0, psi_deriv_1*psi_deriv_1];
    
    % Multiply by Gauss weight
    psi_deriv = psi_deriv * gq.gsw(i);
    
    % Adding onto the matrix
    diffusion_term = diffusion_term + psi_deriv;
end

% Multiply by the Jacobian of the element transformation
diffusion_term = diffusion_term * (D/J);

% diffusion_term = LaplaceElemMatrix(D, eID, msh);

%% Solving reaction term
% Initialise necessary variables
N = 2;                      % Define the GQ level
reaction_term = zeros(2,2); % Initialise an empty 2x2 LEM matrix

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
    reaction_term = reaction_term + psi;
end
reaction_term = reaction_term * msh.elem(eID).J * lambda;

%% Combine two terms
Kelem = diffusion_term - reaction_term;