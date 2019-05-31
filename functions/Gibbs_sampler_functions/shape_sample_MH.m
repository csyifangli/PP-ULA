function [shape_new,accept_shape,wind,walk_shape] = shape_sample_MH(shape_old,scale,x,walk_shape,wind,ii,n0,z,k)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
% 
% This code has been provided by Ningning Zhao (2014/06/06) and commented by
% Marie-Caroline Corbineau in 05/2019.
%
% Input: shape_old: old GGD shape parameter of region k
%        scale: GGD scale parameter of region k
%        x: TRF sample
%        walk_shape: random walk size
%        wind: acceptance window
%        ii: iteration number
%        n0: number of burn-in iterations
%        z: segmentation labels
%        k: region number
%
% Output: shape_new: new shape parameter
%         accept_shape: MH acceptance rate
%         wind: acceptance window
%         walk_shape: random walk size
%
% This function samples the shape parameter using a Metropolis Hastings
% random walk method.
%====================================================================

theta   = x(z==k); 
N       = numel(theta);
shape_cand = shape_old+walk_shape*randn;

% negative log of the posterior distribution w.r.t. shape, and the
% likelihood function wrt shape is ggd
p_func = @(sh) N*log(2*scale^(1/sh)*gamma(1+1/sh))+norm(theta(:),sh)^sh/scale;
p_cand = p_func(shape_cand);
p_old  = p_func(shape_old);

if shape_cand<=3 && shape_cand>=0.01 
    accept = min(exp(p_old-p_cand),1);
else
    accept = 0;
end
    
urn = rand;
if urn<=accept
    shape_new = shape_cand;
    label     = 1;
else
    shape_new = shape_old;
    label     = 0;
end

% adjust the random walk size 
nwind        = 30;
wind         = [wind(1,2:nwind),label];
accept_shape = mean(wind);
if ii-1<=n0
    if  accept_shape<0.4
        walk_shape = walk_shape*(1-0.2);
    elseif accept_shape>0.8
        walk_shape = walk_shape*(1+0.2);
    end
elseif ii-1>n0
    if  accept_shape<0.4
        walk_shape = walk_shape*(1-0.1*exp(-0.01*(ii-n0)));
        elseif accept_shape>0.8
        walk_shape = walk_shape*(1+0.1*exp(-0.01*(ii-n0)));
    end
end
end




