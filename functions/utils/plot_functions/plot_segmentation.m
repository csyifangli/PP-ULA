function [] = plot_segmentation(list,refx,radii,centers,ang)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  list: results, list of cells, each cell is a structure
%         refx: 0: ground-truth is not known, real data
%               1: TRF and segmentation ground-truths are known
%               2: only TRF ground-truth is known
%         radii: radius of each ellipse 
%         centers: center of each ellipse
%         ang: orientation angle of each ellipse
%
% This function plot segmentation visual result where important locations
% are highlighted with green ellipses.
%====================================================================

nb_subplots = length(list);
if refx==1
    size_fig    = [500 500 nb_subplots*300 200];
    figure('pos',size_fig)
    subplot(1,nb_subplots,1)
    imagesc(list{1}.seg),colormap gray,axis off
    title(list{1}.name)
    j = 0;
else
    j = -1;
    size_fig    = [500 500 (nb_subplots-1)*300 200];
    figure('pos',size_fig)
end
for i=2:length(list)
    subplot(1,nb_subplots+j,i+j)
    imagesc(list{i}.seg),colormap gray,axis off
    ellipse(radii(:,1),radii(:,2),ang,centers(:,1),centers(:,2),'g')
    title(list{i}.name)
end 
end

