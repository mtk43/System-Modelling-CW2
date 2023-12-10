function [LEM] = LEMgeneratorLinear(gq, N, basis_type)
%%Function to generate the local element matrix from the

LEM = zeros(2);
% Loop through Gauss points
for it = 1:N
    % Evaluate basis function, either for its value or the gradient
    switch basis_type
        case 'Value'            
            psi_0 = EvalBasis(0, gq.xipts(it));
            psi_1 = EvalBasis(1, gq.xipts(it));
        case 'Gradient'
            psi_0 = EvalBasisGrad(0, gq.xipts(it));
            psi_1 = EvalBasisGrad(1, gq.xipts(it));
    end
    
    % Generate matrix with psi values in correct indices
    psi = [psi_0*psi_0, psi_0*psi_1; psi_1*psi_0, psi_1*psi_1];
    
    % Multiply by Gauss weight
    psi = psi * gq.gsw(it);
    
    % Adding onto the LEM
    LEM = LEM + psi;
end
end