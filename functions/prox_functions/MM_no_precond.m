function[px,bool] = MM_no_precond(x,shape,alpha,precision,label,K)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  x: point at which the proximity operator must be computed
%         shape: GGD shape parameters
%         alpha: GGD scale parameters
%         precision: precision for the iterative solver
%         label: segmentation labels
%         K: number of regions
%
% Output: px: proximity operator
%         bool: 0: MM algorithm did not converge, 1: converged
%
% This function computes the following proximity operator using a
% Majorization-Minimization method:
% argmin_z   1/2||x-z||^2 + sum_k alpha(k)|z_k|^shape_{label_k}
%====================================================================

bool = 1;    % boolean to check if the MM solver has converged
imax = 300;  % maximal number of iterations

% initialization
px      = x;
pxold   = x;
coeff   = {};
objref = CalculObj(px);

for i=1:imax
    %%% MM principle
    for k=1:K
        if shape(k)<1
            pxk                     = px(label{k});
            coeffk                  = zeros(length(label{k}),1);
            coeffk(abs(pxk)>=1e-11) = shape(k).*exp((shape(k)-1).*log(abs(pxk(abs(pxk)>=1e-11))));
            coeff{k}                = coeffk;
        end
    end
    %%% convex proximity operator |.|^p for p>=1
    for k=1:K
        if shape(k)>=1
            px(label{k}) = prox_lp(x(label{k}),shape(k),alpha(k));
        else
            px(label{k}) = prox_lp(x(label{k}),1,alpha(k)*coeff{k});
        end
    end
    list_norm(i) = norm(px(:)-x(:));
    list_obj(i)  = CalculObj(px);
    if norm(px-pxold)/norm(pxold)<precision && CalculObj(px)<=objref,break,end
    if bool==0, break,end
    if i==imax
       disp('----------------------------did not converge outloop--------------------------------')
       fprintf('norm %.4e (obj-obj_ref)/objref = %.4e\n',norm(px-pxold)/norm(pxold),(CalculObj(px)-objref)/obj1ref)
       bool = 0;
    end
    pxold = px;
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% objective function which is supposed to be minimized by the proximity
% operator
function[obj] = CalculObj(px)
    obj = 0.5*sum((x(:)-px(:)).^2);
    for kk=1:K
        obj = obj+alpha(kk)*sum(exp(shape(kk).*log(abs(px(label{kk})))));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end