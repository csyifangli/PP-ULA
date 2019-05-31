function [ psnr ] = PeakSignalToNoiseRatio(x,x_true)
%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  x: estimated image
%         x_true: ground-truth
%
% Output: psnr: peak signal-to-noise ratio
%         
% This function computes the peak signal-to-noise ratio (PSNR) between a
% ground-truth image and an estimated one.
%====================================================================
psnr = 10*log10(numel(x(:))*max(max(x_true(:),x(:)))^2/norm(x_true(:)-x(:))^2);
end


