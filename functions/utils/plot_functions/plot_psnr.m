function [] = plot_psnr(list,refx)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  list: list of cells, each cell is a structure
%         refx: 0: ground-truth is not known, real data
%               1: TRF and segmentation ground-truths are known
%               2: only TRF ground-truth is known
%
% This function plots the PSNR with regards to time.
%====================================================================

if refx>0
    size_fig = [500 500 750 240];
    figure('rend','painters','pos',size_fig)
    legend_names = {[]};
    hold on
    for j=2:length(list)
        if isempty(list{j})==0
            plot(list{j}.time,list{j}.psnr(:,1),list{j}.style,'Linewidth',2)
            legend_names{j-1} = list{j}.name;
        end
    end
    for j=2:length(list)
        if isempty(list{j})==0
            plot(list{j}.time(list{j}.n0+1:end),...
                list{j}.psnr(list{j}.n0+1:end,2),...
                strcat(':',list{j}.style),'Linewidth',2)
        end
    end
    hold off
    set(gca,'fontsize',14)
    legend(legend_names,'Interpreter','latex','Fontsize',14)
    ylabel('PSNR (dB)','Interpreter','latex','Fontsize',14)
    xlabel('Time (s)','Interpreter','latex','Fontsize',14)
    set(gca,'TickLabelInterpreter','latex')
end
end

