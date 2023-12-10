function [ psi ] = EvalBasisQuad(lnid,xipt)
%EvalBasis returns value of quadratic Lagrange basis function
%   The funtion takes in the local node ID (lnid = 0 or 1), and the Gauss point
%   (xipt), and then uses the linear Lagrange basis function equation to
%   determine the value of psi, using the lnid value for the sign of basis
%   function  
    switch lnid
        case 0
            psi = (xipt*(xipt-1))/2;
        case 1
            psi = 1 - (xipt^2);
        case 2
            psi = (xipt*(xipt+1))/2;
    end
end

