function [val] = EvalFieldLambda(msh, eID, xipt)
%EvalField determines the value of a field at a specific node
% Takes in the 1D mesh (msh), the field in question (field), the element ID
% (eID) and the Gauss point (xipt) to find the value of the field

% Extract the correspinding global node IDs
GN = msh.elem(eID).nq;
% GN = msh.elem(eID).n;

% Determine the value of beta based on the eID (based on mesh dimensions
% (0,0.01, 12)
if(eID >=1 && eID <= 2)
    beta = 0;
else
    beta = 0.001;
end

% Define gamma value
gamma = 0.02;

% Combine beta and gamma to form final lambda value
lambda_val = -beta - gamma;

val = 0;
% Find field value for each and sum
for i = 1:length(GN)
    psi = EvalBasisQuad(i-1, xipt); %Evaluate basis function at specified Gauss point
    Lambda = psi * lambda_val;     %Determine local field value
    val = val + Lambda;              %Sum
end


end