function [Melem] = LocalMassMatrix(eID, msh, gq, type)
%Function which generates the local element mass matrix (Melem), of the linear
%reaction term
%%Takes in the element ID (eID), the 1D mesh (msh),and the Gaussian Quadrature scheme (GQ)
%% Initialise necessary variables
N = gq.npts;              % Define the GQ level

%% Evaluate basis functions - linear or quadratic
switch type
    case 'Linear'
        LEM = LEMgeneratorLinear(gq, N, 'Value');
    case 'Quadratic'
        LEM = LEMgeneratorQuad(gq, N, 'Value');        
end

Melem = LEM * msh.elem(eID).J;
end