function [] = plot_TRF_Bmode(list,refx,cst,ulc1,ulc2,width)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  list: results, list of cells, each cell is a structure
%         refx: 0: ground-truth is not known, real data
%               1: TRF and segmentation ground-truths are known
%               2: only TRF ground-truth is known
%         cst: contrast used to compute B-mode images
%         ulc1,ulc2: position of upper-left-hand corner of the two windows
%                    which are used to compute the CNR
%         width: width of the two windows
%
% This function plots the B-mode of the input TRF and RF images. 
% For the RF image this function also plots the windows that are used to 
% compute the contrast-to-noise ratio.
%====================================================================


nb_subplots = length(list);
j           = 0;
if refx>0,j=1; end
size_fig    = [500 500 (nb_subplots+j)*300 200];

figure('pos',size_fig)
%%% plot RF image
subplot(1,nb_subplots+j,1)
imagesc(rf2bmode(list{1}.RF,cst)),colormap gray; axis off    
hold on
rectangle('Position',[ulc1(2) ulc1(1) width(2) width(1)],'Linewidth',1,'Edgecolor','b')
rectangle('Position',[ulc2(2) ulc2(1) width(2) width(1)],'Linewidth',1,'Edgecolor','b')
hold off
title('RF image')
if refx>0
    %%% plot ground-truth TRF
    subplot(1,nb_subplots+j,2)
    imagesc(rf2bmode(list{1}.TRF,cst)),colormap gray; axis off
    title(list{1}.name)
end
for i=2:length(list)
    %%% plot TRF images
    subplot(1,nb_subplots+j,i+j)
    imagesc(rf2bmode(list{i}.TRF,cst)),colormap gray; axis off
    title(list{i}.name)
end
end

