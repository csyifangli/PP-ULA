function [] = print_metrics(image,list,ulc1,ulc2,width,cst,K,refx)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  image: simulation name
%         list: list of cells with results, each cell is a structure
%         ulc1,ulc2: position of upper-left-hand corner of the two windows
%                    which are used to compute the CNR
%         width: width of the two windows
%         cst: constrast used to show the B-mode images
%         K: number of regions
%         refx: 0: ground-truth is not known, real data
%               1: TRF and segmentation ground-truths are known
%               2: only TRF ground-truth is known
%         
% This function prints evaluation metrics (PSNR, SSIM, OA...).
%====================================================================

    fprintf('======== %s ========\n',image)
    if refx>0
        fprintf('-------- PSNR ---- SSIM ---- CNR --------\n')
        for j=2:length(list)
            if isempty(list)==0
                fprintf(strcat(list{j}.name,'   %.1f      %.2f      %.2f\n'),...
                    PeakSignalToNoiseRatio(list{j}.TRF,list{1}.TRF),...
                    ssim(list{j}.TRF,list{1}.TRF),...
                    ContrastToNoiseRatio(rf2bmode(list{j}.TRF,cst),ulc1,ulc2,width));
            end
        end
    else
        fprintf('-------- CNR --------\n')
        for j=2:length(list)
            if isempty(list)==0
                fprintf(strcat(list{j}.name,'   %.2f\n'),...
                    ContrastToNoiseRatio(rf2bmode(list{j}.TRF,cst),ulc1,ulc2,width));
            end
        end
    end
    
    if refx==1
        fprintf('-------- OA --------\n')
        for j=2:length(list)
            if isempty(list)==0
                fprintf(strcat(list{j}.name,'    %.1f\n'),...
                    100*OverallAccuracy(list{j}.seg,list{1}.seg));
            end
        end
    end
    
    fprintf('----------------------\n')
    fprintf('Size %d x %d\n',size(list{1}.RF,1),size(list{1}.RF,2))
    
    fprintf('-------- Time --------\n')
    for j=2:length(list)
        if isempty(list)==0
            fprintf(strcat(list{j}.name,'        %.0f min\n'),list{j}.time(end)/60);
        end
    end
    
    fprintf('-------- MSJ --------\n')  
    for j=2:length(list)
        if isempty(list)==0
            fprintf(strcat(list{j}.name,'     %.1f\n'),list{j}.MSJ(1));
        end
    end

    fprintf('-------- MSJ/s --------\n')
    for j=2:length(list)
        if isempty(list)==0
            fprintf(strcat(list{j}.name,'     %.1f\n'),list{j}.MSJ(2));
        end
    end

    fprintf('================================== GGD parameters ===================================\n') 
    fprintf('         sig2       alpha_1       beta_1       alpha_2       beta_2')
    for ii=3:K
        fprintf('       alpha_%d       beta_%d',ii,ii)
    end
    if refx==1
        %%% True
        fprintf('\nTrue     %.3f      %.3f         %.3f         %.3f         %.3f',list{1}.sig2,...
            list{1}.shape(1),list{1}.scale(1),list{1}.shape(2),list{1}.scale(2))
        for ii=3:K
            fprintf('         %.3f         %.3f',list{1}.shape(ii),list{1}.scale(ii))
        end
    end
    %%% PP-ULA and P-ULA
    for j=2:length(list)
        if isempty(list)==0
            fprintf(strcat('\n',list{j}.name,'   %.3f      %.3f         %.3f       %.3f       %.3f'),...
                list{j}.sig2,list{j}.shape(1),list{j}.scale(1),list{j}.shape(2),list{j}.scale(2));
            for ii=3:K
                fprintf('         %.3f         %.3f',list{j}.shape(ii),list{j}.scale(ii))
            end
        end
    end
    fprintf('\n==================================================================================\n')