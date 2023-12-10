function [val] = EvalFieldDiff(msh, eID, xipt, param)
%EvalField determines the value of a field at a specific node
% Takes in the 1D mesh (msh), the field in question (field), the element ID
% (eID) and the Gauss point (xipt) to find the value of the field

% Extract the correspinding global node IDs
GN = msh.elem(eID).nq;
% GN = msh.elem(eID).n;

% Determine the value of D based on the eID (based on mesh dimensions
% (0,0.01, 12)
% if(eID >=1 && eID <= 2)
%     d = 4e-6;
% elseif(eID >= 3 && eID <= 6)
%     d = 5e-6;
% else
%     d = 2e-6;
% end

if (eID >= 1 && eID <= param.layers.E)
    % Epidermis
    d = param.material.D(1);
elseif (eID > param.layers.E && eID <= param.layers.D)
    % Dermis
    d = param.material.D(2);
else 
    % Sub-cutaneous
    d = param.material.D(3);
end
    

val = 0;
% Find field value for each and sum
for i = 1:length(GN)
    psi = EvalBasisQuad(i-1, xipt); %Evaluate basis function at specified Gauss point
    D = psi * d;                    %Determine local field value
    val = val + D;                  %Sum
end


end