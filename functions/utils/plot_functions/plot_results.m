function [] = plot_results(image,path_results)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  image: simulation name
%         path_results: path to the folder containing simulation results
%
% This function plots results corresponding to a simulation (segmentation, 
% B-mode of the TRF...)
%====================================================================

folder       = strcat(path_results,'/',image);
listing_full = dir(folder);
for ii=3:length(listing_full)
    listing{ii-2} = listing_full(ii).name;
end

%%% load data
noise_seed = 3;
rng(noise_seed);
[y, rescale, remove_padding, refx, K, ~, ~, struct] = input_images(image);
y   = remove_padding(rescale(y));
cst = struct.cst;
if refx>0
    true_TRF = struct.TRF;
    if refx==1
        true_sig2  = struct.sig2;
        true_shape = struct.shape;
        true_scale = struct.scale;
        true_seg   = struct.mask;
        % rearrange segmentation so that comparison with groundtruth is correct
        switch image
            case 'Simu1'
                l = 21;
            case 'Simu2'
                l = 321;
        end
        [true_scale,true_shape] = permutation_GGD(true_scale,true_shape,K,l);
        true_seg = permutation_segmentation(true_seg,K,l);
    end
end

%%% load positions of the windows used to compute the CNR, and the
% characteristics of the ellipses to plot on the segmentation
[ulc1,ulc2,width,centers,radii,ang] = visual_elements(image);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Data RF image
ii=1;
list{ii}.RF  = y;
if refx>0,  end
%%% Ground-truth
if refx>0
    list{ii}.RF_wo_noise = remove_padding(struct.RF_wo_noise);
    list{ii}.name        = 'True';
    list{ii}.TRF   = true_TRF;
    if refx==1    
        list{ii}.style = 'g';
        list{ii}.seg   = true_seg;
        list{ii}.shape = true_shape;
        list{ii}.scale = true_scale;
        list{ii}.sig2  = true_sig2;
    end
end

%%% PP-ULA    
ii=2;
simu_name = startsWith(listing,'Precond1');
if sum(simu_name)==0
    list{ii} = [];
    disp('There are no result corresponding to PP-ULA in the indicated folder.')
elseif sum(simu_name)>1
    list{ii} = [];
    disp('There are several results corresponding to PP-ULA in the indicated folder.')
    disp('Please remove unnecessary results.')
else
    simu_name = listing{simu_name};
    load(strcat(folder,'/',simu_name),'results');
    list{ii}.name  = 'PP-ULA';
    list{ii}.style = 'b';
    list{ii}.TRF   = results.TRF;
    list{ii}.seg   = results.seg;
    list{ii}.shape = results.shape.mmse;
    list{ii}.scale = results.scale.mmse;
    list{ii}.sig2  = results.sig2.mmse;
    list{ii}.time  = results.time_vec;
    list{ii}.MSJ   = results.MSJ;
    if refx>0
        list{ii}.psnr  = results.psnr;
    end
    list{ii}.n0   = results.n0;
    list{ii}.TRF_neglog   = results.TRF_neglog;
    clear('results')
end

%%% P-ULA 
ii=3;
simu_name = startsWith(listing,'Precond0');
if sum(simu_name)==0
    disp('There are no result corresponding to P-ULA in the indicated folder.')
elseif sum(simu_name)>1
    disp('There are several results corresponding to P-ULA in the indicated folder.')
    disp('Please remove unnecessary results.')
else
    simu_name = listing{simu_name};
    load(strcat(folder,'/',simu_name),'results');
    list{ii}.name  = 'P-ULA';
    list{ii}.style = 'c';
    list{ii}.TRF   = results.TRF;
    list{ii}.seg   = results.seg;
    list{ii}.shape = results.shape.mmse;
    list{ii}.scale = results.scale.mmse;
    list{ii}.sig2  = results.sig2.mmse;
    list{ii}.time  = results.time_vec;
    list{ii}.MSJ   = results.MSJ;
    if refx>0
        list{ii}.psnr  = results.psnr;
    end
    list{ii}.n0   = results.n0;
    list{ii}.TRF_neglog   = results.TRF_neglog;
    clear('results')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot segmentation
plot_segmentation(list,refx,radii,centers,ang)

%%%% Plot Bmode of TRF
plot_TRF_Bmode(list,refx,cst,ulc1,ulc2,width)

%%%% Plot distributions
plot_distributions(list,K,refx)

%%%% Print metrics
print_metrics(image,list,ulc1,ulc2,width,cst,K,refx)

%%%% Plot PSNR
plot_psnr(list,refx)

%%%% Plot conditional proba of TRF
plot_conditional_proba_TRF(list,refx,K)
end

