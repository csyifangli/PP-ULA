function modeB = rf2bmode(RF, increase)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  RF: radio-frequency ultrasound image
%         increase: contrast
%
% Output: modeB: B-mode of the RF image
%
% This function converts a RF image into a mode B one using the Hilbert
% transform.
%====================================================================

clear modeB

%%% Initialization
if nargin==1
    increase=0;
elseif nargin<1 || nargin>2
    error('There must be 1 or 2 arguments.')
end

%%% Computation
for i=1:size(RF,3)
    modeB_temp = 20*log10(abs(hilbert(RF(:,:,i)))+increase);
    modeB_temp = modeB_temp-min(modeB_temp(:));
    max_modeB  = max(max(modeB_temp));
    modeB_temp = modeB_temp/max_modeB;
    modeB(:,:,i)=modeB_temp;
end