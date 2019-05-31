function [seg] = permutation_segmentation(seg,K,user_input)
%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  seg: segmentation 
%         K: number of regions
%         user_input: permutation order
%                     For instance, if there are three regions and if 
%                     user_input = 321 then the pixels associated with 
%                     segmentation labels 1, 2 and 3 are now associated 
%                     with labels 3, 2 and 1, respectively.
%
% Output: scale,shape: permuted scale and shape parameters
%
% This function performs a permutation of the segmentation labels when the 
% labels need to be changed (for comparison, for instance).
%====================================================================

user_input = num2str(user_input);
for ii=1:K
    new_ind(ii)  = str2double(user_input(ii));
    seg(seg==ii) = ii*100;
end
for ii=1:K
    seg(seg==ii*100) = new_ind(ii);
end
end