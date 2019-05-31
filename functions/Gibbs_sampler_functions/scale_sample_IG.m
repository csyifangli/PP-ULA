function [scale_new] = scale_sample_IG(x,shape,scale,z,k)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
% 
% This code has been provided by Ningning Zhao and modified by
% Marie-Caroline Corbineau in 05/2019.
%
% Intput: x: TRF sample
%         shape: GGD shape parameter of region k
%         scale: GGD scale parameter of region k
%         z: segmentation labels
%         k: region number
%
% Output: scale_new: new scale parameter of region k
%
% This function samples the scale parameters according to the inverse Gamma
% distribution.
%====================================================================

theta = x(z==k);
N     = numel(theta);

alpha = N/shape;
beta  = norm(theta(:),shape)^shape;

scale_new = 1/gamrnd(alpha,1/beta); 
if scale_new<=0 
    scale_new = scale;
end

end