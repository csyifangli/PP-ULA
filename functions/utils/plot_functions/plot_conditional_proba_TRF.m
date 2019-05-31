function [] = plot_conditional_proba_TRF(list,refx,K)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  list: results, list of cells, each cell is a structure
%         refx: 0: ground-truth is not known, real data
%               1: TRF and segmentation ground-truths are known
%               2: only TRF ground-truth is known
%         K: number of regions
%
% This function plots the negative log of the conditional distribution
% of the TRF (Equation (10) in the paper):
% 1/(2*sigma^2)||y-Hx||^2 + sum_k |x_k|^{shape_k}/scale_k
%====================================================================

size_fig = [500 500 700 240];
figure('pos',size_fig)

if refx==1
    % compute the negative log of the conditional distribution
    % of the true TRF 
    true_TRF_neglog = 0;
    for k=1:K
        true_TRF_neglog = true_TRF_neglog + sum(abs(list{1}.TRF(list{1}.seg==k)).^list{1}.shape(k))/list{1}.scale(k);
    end
    true_TRF_neglog = true_TRF_neglog + 0.5*norm(list{1}.RF-list{1}.RF_wo_noise,'fro')^2/list{1}.sig2;
end
    
hold on
for j=2:length(list)
    plot(list{j}.time,list{j}.TRF_neglog,list{j}.style,'Linewidth',2)
    legend_names{j-1} = list{j}.name;
end
if refx == 1 
    if length(list)==2
        len = ceil(list{2}.time(end));
    else
        len = ceil(list{3}.time(end));
    end
    plot(1:len,true_TRF_neglog.*ones(1,len),'g','Linewidth',2)
    legend_names{length(legend_names)+1} = 'Ground-truth';
end
hold off
set(gca,'fontsize',14)
xlabel('Time (s)','Interpreter','latex','Fontsize',14)
set(gca,'TickLabelInterpreter','latex')
legend(legend_names,'Interpreter','latex','Fontsize',14)
title('$-\log(p(x|y,\sigma^2,\alpha,\beta,z))$','Interpreter','latex','Fontsize',15)
end

