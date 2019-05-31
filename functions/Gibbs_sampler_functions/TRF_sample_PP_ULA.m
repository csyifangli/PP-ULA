function [x_sample,VV,bool] = TRF_sample_PP_ULA(VV,z,scale,shape,sig2,x_sample,...
    alpha,precision,Q_1HtH,Q_1,Q_1Hty,Q_12,norm_Q_1,m,n,K,precond_bool,HtH,Hty)
        
%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  VV: dual variable used for warm-restart
%         z: segmentation labels
%         scale: GGD scale parameters
%         shape: GGD shape parameters
%         sig2: noise variance
%         x_sample: previous sample of the TRF
%         alpha: stepsize is alpha*sig2 (we use alpha by simplification
%         sig2/sig2)
%         precision: precision for Majorization-Minimization subiterations
%         Q_1HtH: Q^{-1}*H^T*H with Q the variable metric
%         Q_1: inverse of the variable metric
%         Q_1Hty: Q^{-1}*H^T*y
%         Q_12: Q^{1/2}
%         norm_Q_1: norm of the inverse of the variable metric
%         m: number of lines
%         n: number of columns
%         K: number of regions (labels)
%         precond_bool: 0: without variable metric (P-ULA), 1: with
%                       variable metric (PP-ULA)
%         HtH: H^T*H
%         Hty: H^T*y
%
% Output: x_sample: next TRF sample
%         VV: dual variable used for warm-restart
%         bool: 0: MM algorithm did not converge, 1: converged
%
% This function outputs the next sample of the TRF using PP-ULA (forward-
% backward + variable metric) or P-ULA (forward-backward). Equation (12) 
% in the paper.
%====================================================================

for kk=1:K
    label{kk} = find(z==kk);
end
switch precond_bool
    case 1
        %preconditioned forward
        xtemp            = x_sample-alpha.*(real(ifft2(Q_1HtH.*fft2(x_sample)))-Q_1Hty);
        %preconditioned backward
        [x_mean,VV,bool] = prox_lp_precond(xtemp(:),shape,sig2*alpha./scale,precision,Q_1,label,norm_Q_1,VV,m,n,K);
        perturbation     = sqrt(2*alpha*sig2).*real(ifft2(Q_12.*fft2(normrnd(zeros(m,n),ones(m,n)))));
    case 0
        %forward
        xtemp         = x_sample-alpha.*(real(ifft2(HtH.*fft2(x_sample)))-Hty);
        %backward
        [x_mean,bool] = MM_no_precond(xtemp(:),shape,sig2*alpha./scale,precision,label,K);
        perturbation  = sqrt(2*alpha*sig2).*normrnd(zeros(m,n),ones(m,n));
        VV            = [];
end
%Euler discretization of Langevin diffusion
x_sample = reshape(x_mean',m,n) + perturbation;
end

