function [percent] = compute_isolated_points(seg)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  seg: segmented image
%
% Output: percent: proportion of isolated points
%
% This function computes the proportion of isolated points in the 
% segmentation using a 3 × 3 median filter.
%====================================================================

seg_filt          = medfilt2(seg);
seg_filt(1,1)     = floor(median([seg(1:2,1);seg(1:2,2)]));
seg_filt(1,end)   = floor(median([seg(1:2,end-1);seg(1:2,end)]));
seg_filt(end,1)   = floor(median([seg(end-1:end,1);seg(end-1:end,2)]));
seg_filt(end,end) = floor(median([seg(end-1:end,end-1);seg(end-1:end,end)]));
percent           = sum(sum(abs(seg_filt-seg)>0))/numel(seg);
end
