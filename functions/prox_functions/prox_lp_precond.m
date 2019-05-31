function[px,v,bool] = prox_lp_precond(x,shape,alpha,precision,A_1,label,norm_A_1,v,m,n,K)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  x: point at which the proximity operator must be computed
%         shape: GGD shape parameters
%         alpha: multiplicative factor
%         precision: precision for the DFB subiterations
%         A_1: inverse of variable metric
%         label: segmentation labels
%         norm_A_1: norm of the inverse of variable metric
%         v: dual variable used for warm-restart
%         m: number of lines
%         n: number of columns
%         K: number of regions (labels)
%
% Output: px: proximity operator
%         v: dual variable used for warm-restart
%         bool: 0: algorithm did not converge, 1: converged
%
% This function computes the following proximity operator,
% argmin_z   1/2||x-w||_Q^2 + lam*sum_k(|w_k|^p)
% using a Majorization-Minimization method combines with the dual forward-
% backward algorithm.
%
% Reference for the dual forward-backward algorithm:
% P. L. Combettes, D. Dung, and B. C. Vu, “Proximity for sums of com-
% posite functions,” Journal of Mathematical Analysis and applications,
% vol. 380, no. 2, pp. 680–688, 2011
%====================================================================

bool = 1;    % boolean to check convergence
imax = 300;  % max number of iteration for forward backward
jmax = 300;  % max number of iterations for dual FB (prox computation)
    
if (isempty(v))
   v = zeros(size(x)); % dual variable
end
rho   =  1/(norm_A_1);
delta =  min(1,rho) - 1e-8;
step  =  1.99*rho-delta ; % should be less than 2*rho-delta

%initialization
px      = x;
xold    = x;
coeff   = {};
coeffb  = {};
obj1ref = CalculObj1(px);

for i=1:imax
    %%% MM principle
    for k=1:K
        if shape(k)<1
        pxk                     = px(label{k});
        coeffk                  = zeros(length(label{k}),1);
        coeffk(abs(pxk)>=1e-11) = shape(k).*exp((shape(k)-1).*log(abs(pxk(abs(pxk)>=1e-11))));
        coeff{k}                = coeffk;
        coeffbk                 = zeros(length(label{k}),1);
        coeffbk(abs(pxk)>=1e-11)= (1-shape(k)).*exp(shape(k).*log(abs(pxk(abs(pxk)>=1e-11))));
        coeffb{k}               = coeffbk;
        end
    end
    %%% start DualFB algorithm 
    pxold  = px;
    objref = CalculObj(px,coeff,coeffb);
    for j=1:jmax
        px = x - reshape(real(ifft2(A_1.*fft2(reshape(v,m,n)))),m*n,1);
        u  = v + step *px;  
        for k=1:K
            if shape(k)>=1
               v(label{k})  =  u(label{k}) - step .*prox_lp(u(label{k})./step,shape(k),alpha(k)./step);
            else
               v(label{k})  =  u(label{k}) - step .*prox_lp(u(label{k})./step,1,alpha(k)*coeff{k}./step);
            end
        end
        
        %%% sanity check about stopping criterion of DualFB
        if j>1
           if norm(px-pxold)/norm(pxold)<precision && CalculObj(px,coeff,coeffb)<=objref,break,end
           if j==jmax
                if norm(px-pxold)/norm(pxold)<precision && CalculObj(px,coeff,coeffb)>objref
                    fprintf('starting point is best because of objective (obj-obj_ref)/objref=%.4e \n',(CalculObj(px,coeff,coeffb)-objref)/objref)
                    px=xold;
                    break
                else
                    disp('did not converge because of precision\n')
                end
                bool = 0;
           end
        end
        pxold = px;
    end
    %%% sanity check about stopping criterion of MM method
    if norm(px-xold)/norm(xold)<precision && CalculObj1(px)<=obj1ref,break,end
    if bool==0,break,end
    if i==imax
       disp('----------------------------did not converge outloop--------------------------------')
       fprintf('norm %d \n',norm(px-xold)/norm(xold))
       bool = 0;
    end
    xold = px;
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% objective function to be minimized by DualFB
function[obj] = CalculObj(px,coeff,coeffb)
    obj = 0.5*sum((x-px).*reshape(real(ifft2(fft2(reshape(x-px,m,n))./A_1)),m*n,1));
    for kk=1:K
        if shape(kk)>=1
            obj = obj+alpha(kk)*sum(exp(shape(kk).*log(abs(px(label{kk})))));
        else
            obj = obj+alpha(kk)*sum(coeff{kk}.*abs(px(label{kk}))+coeffb{kk});
        end
    end
end
% objective function to be minimized by MM algorithm
function[obj] = CalculObj1(px)
    obj = 0.5*sum((x-px).*reshape(real(ifft2(fft2(reshape(x-px,m,n))./A_1)),m*n,1));
    for kk=1:K
        obj = obj+alpha(kk)*sum(exp(shape(kk).*log(abs(px(label{kk})))));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

    

