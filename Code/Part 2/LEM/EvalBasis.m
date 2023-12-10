function [ psi ] = EvalBasis(lnid,xipt)
%EvalBasis returns value of linear Lagrange basis function
%   The funtion takes in the local node ID (lnid = 0 or 1), and the Gauss point
%   (xipt), and then uses the linear Lagrange basis function equation to
%   determine the value of psi, using the lnid value for the sign of basis
%   function

    sign = (-1)^(lnid+1);       %Determine sign of basis function based on local node
    psi = (1 + sign*xipt)/2;    %Calculate basis function value
end