function [CNR] = ContrastToNoiseRatio(Bmode,ulc1,ulc2,width)
%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  x: image
%         ulc1,ulc2: upper left hand corner of windows 1 and 2
%         width:  size of the windows
%         
% Output: CNR: contrast-to-noise ratio 
%         
% This function computes the contrast-to-noise ratio between two windows in
% an image.
% CNR = |\mu_1-\mu_2|/\sqrt(\sigma_1^2+\sigma_2^2)
% mu_1,2: average value in region 1,2
% sigma_1,2: standard deviation in region 1,2
%====================================================================

win1 = Bmode(ulc1(1):ulc1(1)+width(1),ulc1(2):ulc1(2)+width(2));
win2 = Bmode(ulc2(1):ulc2(1)+width(1),ulc2(2):ulc2(2)+width(2));
CNR  = abs(mean(win1(:))-mean(win2(:)))/sqrt(var(win1(:))+var(win2(:)));
end

