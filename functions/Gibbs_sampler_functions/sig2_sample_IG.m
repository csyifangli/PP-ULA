function sig2 = sig2_sample_IG(x,DFT,y)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
% 
% This code has been provided by Ningning Zhao and commented by 
% Marie-Caroline Corbineau 05/2019.
%
% Intput: x: TRF sample
%         DFT: direct Fourier transform of the PSF
%         y: RF image
%
% Output: sig2: sample of the noise variance for current iteration
%
% This function samples the noise variance according to its conditional 
% posterior distribution which is the inverse Gamma function.
% p(sig2|...) \sim IG(N/2,1/2*||y-Hx||^2)
%====================================================================

if nargin < 3
    error('Not enough inputs');
end

temp  = y-real(ifft2(DFT.*fft2(x)));
sig2  = 1/gamrnd(numel(x)/2,1/(0.5*norm(temp(:))^2));    

end