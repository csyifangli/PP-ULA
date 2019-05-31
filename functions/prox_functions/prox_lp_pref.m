function pr = prox_lp_pref(x,p,gamma)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  x: point at which the proximity operator must be computed
%         p: exponent belonging to {1,4/3,3/2,2,3,4}
%         gamma: multiplicative factor
%
% Output: pr: proximity operator
%
% This function computes the following proximity operator when p takes
% its value in {1,4/3,3/2,2,3,4}, for which the solution has a closed form.
% argmin_z   1/2||x-z||^2 + gamma*sum_k(|z_k|^p)
%
% Reference: C. Chaux, P. L. Combettes, J.-C. Pesquet, and V. R. Wajs, 
% “A variational formulation for frame-based inverse problems,” 
% Inverse Problems, vol.23, no. 4, pp. 1495–1518, June 2007.
%====================================================================

 gamma = gamma+eps;
 switch p
 case 1
     pr = max(abs(x)-gamma,0).*sign(x);
 case 4/3
     pr = sqrt((256/729)*gamma.^3+x.^2);
     pr = x+4*gamma./3./(2^(1/3)).*((pr-x).^(1/3)-(pr+x).^(1/3));
 case 3/2
     pr = x+9/8*gamma.*sign(x).*(gamma-sqrt(gamma.^2+16/9*abs(x)));
 case 2
     pr = x./(1+2.*gamma);
 case 3
     pr = sign(x).*(sqrt(1+12*gamma.*abs(x))-1)./(6*gamma);
 case 4
     pr = sqrt(x.^2+1./(27*gamma));
     pr = 1/(8*gamma).^(1/3).*((x+pr).^(1/3)-(-x+pr).^(1/3));
 otherwise
     error('unknown prox');
 end
