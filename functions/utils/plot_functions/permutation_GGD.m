function [scale,shape] = permutation_GGD(scale,shape,K,user_input)
%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  scale,shape: scale and shape parameters
%         K: number of regions
%         user_input: permutation order
%                     For instance, if there are three regions and if 
%                     user_input = 321 then the GGD parameters associated
%                     with segmentation labels 1, 2 and 3 are now
%                     associated with labels 3, 2 and 1, respectively.
%
% Output: scale,shape: permuted scale and shape parameters
%
% This function performs a permutation of the GGD parameters when the
% labels need to be changed (for comparison, for instance).
%====================================================================

user_input  = num2str(user_input);
temp_scale  = scale;
temp_shape  = shape;
for ii=1:K
    new_ind(ii)        = str2num(user_input(ii));
    scale(new_ind(ii)) = temp_scale(ii);
    shape(new_ind(ii)) = temp_shape(ii);
end
end
