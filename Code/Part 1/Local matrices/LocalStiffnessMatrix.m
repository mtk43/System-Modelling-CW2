function [Kelem] = LocalStiffnessMatrix(eID, msh, D, lambda, gq, type)
%%This function uses constant coefficient values (should definitely be
%%generalised later)

% Extract J value
J = msh.elem(eID).J;

% Define gq level
N = gq.npts;

switch type
    case 'Linear'
        %% Linear
        % Diffusion term
        LEM_diffusion = LEMgeneratorLinear(gq, N, 'Gradient');
        diffusion_term = LEM_diffusion * (D/J);

        % Reaction term
        LEM_reaction  = LEMgeneratorLinear(gq, N, 'Value');
        reaction_term = LEM_reaction * msh.elem(eID).J * lambda;
        
    case 'Quadratic'
        %% Quadratic
        % Diffusion term
        LEM_diffusion = LEMgeneratorQuad(gq, N, 'Gradient');
        diffusion_term = LEM_diffusion * (D/J);

        % Reaction term
        LEM_reaction  = LEMgeneratorQuad(gq, N, 'Value');
        reaction_term = LEM_reaction * msh.elem(eID).J * lambda;
end
%% Combine two terms
Kelem = diffusion_term - reaction_term;