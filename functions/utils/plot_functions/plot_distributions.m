function [] = plot_distributions(list,K,refx)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  list: results, list of cells, each cell is a structure
%         K: number of regions
%         refx: 0: ground-truth is not known, real data
%               1: TRF and segmentation ground-truths are known
%               2: only TRF ground-truth is known
%
% This function plots generalized Gaussian distributions.
%====================================================================

size_fig = [500 500 K*450 250];
figure('pos',size_fig)
for i=1:K
    for j=2:length(list)
        sh_list(j-1) = list{j}.shape(i);
        sc_list(j-1) = list{j}.scale(i);
        if j==2,line_style = '-';else,line_style = ':';end
        names_colors_lin{j-1} = {list{j}.name,list{j}.style,line_style};
    end
    if refx==1
        sh_list(length(sh_list)+1) = list{1}.shape(i);
        sc_list(length(sc_list)+1) = list{1}.scale(i);
        names_colors_lin{length(names_colors_lin)+1} = {list{1}.name,list{1}.style,':'};
    end
    subplot(1,K,i)
    plot_distributions_one_region(sh_list,sc_list,names_colors_lin)
    ylabel(strcat('$p(x|\alpha_',num2str(i),',\beta_',num2str(i),')$'),'Interpreter','latex','Fontsize',14)
    xlabel('$x$','Interpreter','latex','Fontsize',14)
    clear('sh_list','sc_list','names_colors_lin')
end

%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_distributions_one_region(sh,sc,names_colors_lin)
    %histogram fit of generalized Gaussian distribution
    %input: sh: GGD shape parameter
    %       sc: GGD scale parameter 
    %       names_colors_lin: name for legend, colors and linestyles for plots

    %%% to get 98% of the distribution
    fun  = @(z) 1/2.*(1+sign(z).*gammainc(abs(z).^sh(1)./sc(1),1/sh(1)))-0.99;
    xmax = fzero(fun,0);
    x    = -xmax:.01:xmax;
    ymax = 0;
    hold on
    for ii=1:length(sh)
        pdf     = sh(ii)/(2*sc(ii)^(1/sh(ii))*gamma(1/sh(ii)))*exp(-abs(x).^sh(ii)./sc(ii));
        plot(x,pdf,names_colors_lin{ii}{2},'linewidth',2,'linestyle',names_colors_lin{ii}{3})
        leg{ii} = names_colors_lin{ii}{1};
        ymax    = max(max(pdf(:)),ymax);
    end
    hold off
    legend(leg,'Interpreter','latex','Fontsize',12)
    set(gca,'fontsize',13)
    set(gca,'TickLabelInterpreter','latex'), axis tight
    xlim([-xmax xmax])
    ylim([0 1.1*ymax])
end 
%%%%%%%%%%%%%%%%%%%%%%%
end

