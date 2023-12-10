function [Kelem] = LocalStiffnessMatrix(eID, msh, gq, type, param)
%%This function uses constant coefficient values (should definitely be
%%generalised later)

% Extract J value
J = msh.elem(eID).J;

% Define gq level
N = gq.npts;

switch type
    case 'Linear'
%         %% Linear
%         % Diffusion term
%         LEM_diffusion = LEMgeneratorLinear(gq, N, 'Gradient');
%         diffusion_term = LEM_diffusion * (D/J);
% 
%         % Reaction term
%         LEM_reaction  = LEMgeneratorLinear(gq, N, 'Value');
%         reaction_term = LEM_reaction * msh.elem(eID).J * lambda_val;
        
    case 'Quadratic'
        %% Quadratic
        % Diffusion term
%         LEM_diffusion = LEMgeneratorQuad(gq, N, 'Gradient');
          LEM_diffusion = zeros(3);
        
        for i = 1:N
           % Evaluate the reaction coefficient at the current Gauss point
           D_val = EvalFieldDiff(msh, eID, gq.xipts(i), param); 
           
           % Evaluate basis function, either for its value or the gradient
           psi_0 = EvalBasisGradQuad(0, gq.xipts(i));
           psi_1 = EvalBasisGradQuad(1, gq.xipts(i));
           psi_2 = EvalBasisGradQuad(2, gq.xipts(i));
           
           % Generate matrix with psi values in correct indices
           psi = [psi_0*psi_0, psi_1*psi_0, psi_2*psi_0; psi_0*psi_1, psi_1*psi_1, psi_2*psi_1; psi_0*psi_2, psi_1*psi_2, psi_2*psi_2];
           
           % Add weighted contribution of current Gauss point to local
           % matrix
           LEM_diffusion = LEM_diffusion + (D_val * psi * gq.gsw(i));
           
        end        
        diffusion_term = LEM_diffusion / J;

        % Reaction term
%         LEM_reaction  = LEMgeneratorQuad(gq, N, 'Value');
          LEM_reaction = zeros(3);
          
        for i = 1:N
           % Evaluate the reaction coefficient at the current Gauss point
           lambda_val = EvalFieldLambda(msh, eID, gq.xipts(i)); 
           
           % Evaluate basis function, either for its value or the gradient
           psi_0 = EvalBasisQuad(0, gq.xipts(i));
           psi_1 = EvalBasisQuad(1, gq.xipts(i));
           psi_2 = EvalBasisQuad(2, gq.xipts(i));
           
           % Generate matrix with psi values in correct indices
           psi = [psi_0*psi_0, psi_1*psi_0, psi_2*psi_0; psi_0*psi_1, psi_1*psi_1, psi_2*psi_1; psi_0*psi_2, psi_1*psi_2, psi_2*psi_2];
           
           % Add weighted contribution of current Gauss point to local
           % matrix
           LEM_reaction = LEM_reaction + (lambda_val * psi * gq.gsw(i));
           
        end
        reaction_term = LEM_reaction * msh.elem(eID).J;
        
        
end
%% Combine two terms
Kelem = diffusion_term - reaction_term;