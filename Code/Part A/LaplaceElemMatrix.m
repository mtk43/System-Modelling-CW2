function [elemat] = LaplaceElemMatrix(D, eID, msh)

% Determine the element's Jacobian from the Mesh
J = msh.elem(eID).J;

% Create the 2x2 local element matrix
elemat = [D/(2*J), -D/(2*J); -D/(2*J), D/(2*J)];