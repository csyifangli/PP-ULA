function [ x ] = prox_lp( w,p,lam)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  w: point at which the proximity operator must be computed
%         p: exponent >=1
%         lam: multiplicative factor
%
% Output: x: proximity operator
%
% This function computes the following proximity operator.
% argmin_z   1/2||x-w||^2 + lam*sum_k(|w_k|^p) where p>=1
%====================================================================

prec    = 1e-1;
absw    = abs(w);
maxiter = 50;

pref        = [1 4/3 3/2 2 3];
[minp,imin] = min(abs(p-pref));
if minp==0
    %closed-form solution for some values of the exponent
    x = sign(w).*prox_lp_pref(absw,pref(imin),lam); 
else
    %objective function that should be minimized
    CalculObj = @(z) 0.5*sum((w-z).^2) + lam*sum(abs(z).^p);
    obj_ref   = CalculObj(w);
    %we known that the minimum is between x_min and x_max
    x_max     = absw;
    x_min     = zeros(size(absw));
    x_old     = 0.5.*(x_max+x_min);
    %use dichotomic search, the objective is convex 
    for i=1:maxiter
        x              = 0.5.*(x_max+x_min);
        if i>1 && max(abs(x-x_old))<prec && CalculObj(sign(w).*x)<=obj_ref,break,end
        deriv          = x-absw+lam*p*exp((p-1)*log(x));
        x_max(deriv>0) = x(deriv>0);
        x_min(deriv<0) = x(deriv<0);
        x_old          = x;
    end
    x = sign(w).*x;
end 
end

