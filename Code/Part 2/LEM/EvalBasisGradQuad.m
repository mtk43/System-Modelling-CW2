function [ dpsidxi ] = EvalBasisGradQuad(lnid,xipt)
%EvalBasisGrad Returns gradient of basis functions
%   Returns the gradients of the linear Lagrange basis functions for a
%   specificed local node id (lnid = 0 or 1) and xipt (gradient is a constant
%   value in this case, but for higher order basis functions, would vary
%   with xi

    %Use the node id to generate the sign of the basis gradient - ie.
    %either + or -. when lnid=0, sign is -ve, when lnid=1, sign is +ve.
    sign = (-1)^(lnid+1);
    dpsidxi = 0.5 * sign;
    
    switch lnid
        case 0
            dpsidxi = xipt - (1/2);
        case 1
            dpsidxi = -2*xipt;
        case 2
            dpsidxi = xipt + (1/2);
    end

end

