function [LEM] = LEMgeneratorQuad(gq, N, basis_type)

LEM = zeros(3);
% Loop through Gauss points
for it = 1:N
    % Evaluate basis function, either for its value or the gradient
    switch basis_type
        case 'Value'
            psi_0 = EvalBasisQuad(0, gq.xipts(it));
            psi_1 = EvalBasisQuad(1, gq.xipts(it));
            psi_2 = EvalBasisQuad(2, gq.xipts(it));
        case 'Gradient'
            psi_0 = EvalBasisGradQuad(0, gq.xipts(it));
            psi_1 = EvalBasisGradQuad(1, gq.xipts(it));
            psi_2 = EvalBasisGradQuad(2, gq.xipts(it));
    end
    
    % Generate matrix with psi values in correct indices
    psi = [psi_0*psi_0, psi_0*psi_1, psi_2*psi_0; psi_0*psi_1, psi_1*psi_1, psi_2*psi_1; psi_0*psi_2, psi_1*psi_2, psi_2*psi_2];
    
    % Multiply by Gauss weight
    psi = psi * gq.gsw(it);
    
    % Adding onto the LEM
    LEM = LEM + psi;
end
end