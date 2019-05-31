function [ OA ] = OverallAccuracy( z,z_true )
%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  z: estimated segmentation
%         z_true: ground-truth
%
% Output: OA: overall accuracy
%         
% This function computes the proportion of labels that are correctly
% predicted.
%====================================================================

comp = z(:)-z_true(:);
OA   = sum(comp==0)/numel(z_true(:));
end

