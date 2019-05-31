function [results] = demo(image, n0, Nmc, precond_bool, noise_seed,granularity)
%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  image: simulation name
%         n0: number of burn-in iterations
%         Nmc: number of iterations after burn-in
%         precond_bool: 1 for PP-ULA, 0 for P-ULA (without variable metric)
%         noise_seed: seed of random generator for reproducibility
%
% Output: results: structure aggregating results
%                  results.psnr:         peak signal-to-noise ratio 
%                  results.MSJ:          mean square jump (MSJ) and MSJ per second
%                  results.TRF_neglog:   conditional probability of TRF
%                  results.time_vec:     time after each iteration
%                  results.scale.mmse:   MMSE of scale parameters 
%                  results.scale.sample: scale parameters for each iteration
%                  results.shape.mmse:   MMSE of shape parameters
%                  results.shape.sample: shape parameters for each iteration
%                  results.sig2.mmse:    MMSE of noise variance
%                  results.sig2.sample:  noise variance for each iteration
%                  results.seg:          estimated segmentation
%                  results.TRF:          estimated TRF
%                  results.n0:           number of burn-in iterations
%
% This function runs the hybrid Gibbs sampler with PP-ULA or P-ULA 
% to perform joint segmentation and deconvolution of an ultrasound image.
%====================================================================

rng(noise_seed) 
% load data corresponding to the simulation
[y, rescale, remove_padding, refx, K, DFT, psf_ref, struct] = input_images(image);

if refx>0
    TRF_true = struct.TRF;
end
[m,n]      = size(y);
steps      = Nmc+n0;

if refx>0
    compute_psnr = @(x) PeakSignalToNoiseRatio(rescale(remove_padding(x)),TRF_true);
end

% plot data
switch refx
    case 0
        size_fig = [500 500 300 200];
        figure('pos',size_fig)
        imagesc(rf2bmode(rescale(remove_padding(y)),struct.cst))
        colormap gray, axis off, title('RF image')
        fprintf('Real data, image size = %d x %d.\n',size(remove_padding(y),1),size(remove_padding(y),2))
    case 1
        size_fig = [500 500 3*300 200];
        figure('pos',size_fig)
        subplot(131)
        imagesc(rf2bmode(rescale(remove_padding(y)),struct.cst))
        colormap gray, axis off, title('RF image')
        subplot(132)
        imagesc(rf2bmode(TRF_true,struct.cst))
        colormap gray, axis off, title('True TRF image')
        subplot(133)
        imagesc(struct.mask)
        colormap gray, axis off, title('True segmentation')
        fprintf('Simulated data, image size = %d x %d.\n',size(TRF_true,1),size(TRF_true,2))
    case 2
        size_fig = [500 500 2*300 200];
        figure('pos',size_fig)
        subplot(121)
        imagesc(rf2bmode(rescale(remove_padding(y)),struct.cst))
        colormap gray, axis off, title('RF image')
        subplot(122)
        imagesc(rf2bmode(TRF_true,struct.cst))
        colormap gray, axis off, title('True TRF image')
        fprintf('Tissue-mimicking simulated data, image size = %d x %d.\n',size(TRF_true,1),size(TRF_true,2))
end
pause(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% noise variance
sig2_sample    = zeros(1,steps);
sig2_sample(1) = 0.01; % this has no effect since the noise variance is 
                       % the first parameter to be sampled

%%% TRF
x_sample  = deconvwnr(y,psf_ref/sum(psf_ref(:)),10/var(y(:)));
x_hat     = x_sample;

%%% labels z 
filtered  = medfilt2(rf2bmode(x_sample,1),[7 7]);
thresh    = multithresh(filtered,K-1);
z_sample  = imquantize(filtered,thresh);
label_z   = zeros(m,n,K);

%%% shape parameter
shape_sample      = zeros(K,steps);
shape_sample(:,1) = rand(K,1).*1+0.5;

%%% scale parameter 
scale_sample      = zeros(K,steps);
scale_sample(:,1) = rand(K,1)*199+1; 

%%% Metropolis Hastings parameters for GGD shape params
rand(1,30); % for reproducibility of the resutls presented in the paper
nwind             = 30;
wind_shape        = floor(rand(1,nwind)+0.6);
accept_shape      = zeros(K,steps);
walk_shape        = zeros(K,steps);
walk_shape(:,1)   = 1e-10;

%%% Mean Square Jump
MSJ = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% store useful quantities
HtH = abs(DFT).^2;
Ht  = conj(DFT);
Hty = real(ifft2(Ht.*fft2(y)));

%%% PP-ULA parameters
if precond_bool == 1
    alpha_prox  = 0.09;             %should be less than 1 
else
    alpha_prox  = 1.99/max(HtH(:)); %should be less than 2/max(HtH(:))
end
%for computation of the proximity operator
precision = 5e-2;
VV        = [];
%metric
lambda   = 0.1; % to have a positive definite metric 
Q_1      = 1./(HtH+lambda.*ones(m,n));
Q_12     = sqrt(Q_1);
norm_Q_1 = max(max(Q_1));
Q_1HtH   = HtH.*Q_1;
Q_1Ht    = Ht.*Q_1;
Q_1Hty   = real(ifft2(Q_1Ht.*fft2(y)));

TRF_sample_prox = @(VV,z,s,xi,sig2,x) TRF_sample_PP_ULA(VV,z,s,xi,sig2,x,...
            alpha_prox,precision,Q_1HtH,Q_1,Q_1Hty,Q_12,norm_Q_1,m,n,K,precond_bool,HtH,Hty);

%%%%%%%%%%%%%%%%%%%%%%%%%% Hybrid Gibbs Sampler %%%%%%%%%%%%%%%%%%%%%%%%%%

time_vec(1)       = 0;
for ii=2:steps
    if ii>2 && mod(ii-2,100)==0 
        if refx >0
            fprintf('Iterations performed: %d/%d | PSNR: TRF sample = %.2f  TRF hat = %.2f\n',ii-2,steps-1,psnr_vec(ii-1,1),psnr_vec(ii-1,2))
        else
            fprintf('Iterations performed: %d/%d\n',ii-2,steps-1)
        end
    end
 tic
 
 %%% Sample sig2
 sig2_sample(ii) = sig2_sample_IG(x_sample,DFT,y);
 
 %%% Sample the labels of RF image
 z_temp        = z_sample;
 [z_b,ind_b]   = label_sample_MRF_4neighbours(z_temp,x_sample,K,granularity,shape_sample(:,ii-1),scale_sample(:,ii-1),1);
 z_temp(ind_b) = z_b;    
 [z_w,ind_w]   = label_sample_MRF_4neighbours(z_temp,x_sample,K,granularity,shape_sample(:,ii-1),scale_sample(:,ii-1),2);
 z_temp(ind_w) = z_w;
 z_sample      = z_temp;
        
 %%% Sample shape and scale parameters
 for kk=1:K   
     [shape_sample(kk,ii),accept_shape(kk,ii),wind_shape,walk_shape(kk,ii)]=...
        shape_sample_MH(shape_sample(kk,ii-1),scale_sample(kk,ii-1),x_sample,walk_shape(kk,ii-1),wind_shape,ii,n0,z_sample,kk);
     [scale_sample(kk,ii)] = scale_sample_IG(x_sample,shape_sample(kk,ii),scale_sample(kk,ii-1),z_sample,kk);    
 end
         
 %%% Sample the TRF with PP-ULA  
 x_sample_before    = x_sample;
 [x_sample,VV,bool] = TRF_sample_prox(VV,z_sample,scale_sample(:,ii),shape_sample(:,ii),sig2_sample(ii),x_sample);
 if bool==0,break,end

 %%% Calculate MMSE estimates
 if ii<=n0
     x_hat = x_sample;
 else
     %MMSE estimate of the TRF
     x_hat = (x_hat*(ii-1-n0)+x_sample)/(ii-n0);  
     if ii>=n0+2
         % mean square jump
         MSJ = MSJ + norm(x_sample(:)-x_sample_before(:))^2; 
     end
 end

 %MMSE estimate for the segmentation labels
 for kk=1:K
     if ii>n0 
         ind             = find(z_sample==kk);
         temp            = label_z(:,:,kk);
         temp(ind)       = temp(ind)+1;
         label_z(:,:,kk) = temp;
     end
  end

  time_vec(ii) = time_vec(ii-1)+toc;
  
  %%% Compute PSNR and the conditional proba of the TRF
  if refx>0, psnr_vec(ii,:) = [compute_psnr(x_sample)  compute_psnr(x_hat)]; end
  TRF_neglog_vec(ii) = compute_TRF_neglog(x_sample,z_sample,...
      shape_sample(:,ii),scale_sample(:,ii),sig2_sample(ii));
end
if refx >0
    fprintf('Iterations performed: %d/%d | PSNR: TRF sample = %.2f  TRF hat = %.2f\n',steps-1,steps-1,psnr_vec(end,1),psnr_vec(end,2))
else
    fprintf('Iterations performed: %d/%d\n',steps-1,steps-1)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculate the estimators
[~,ind] = max(label_z,[],3);
z_map   = ind;

sig2_est.mmse   = mean(sig2_sample(n0+1:end));
sig2_est.sample = sig2_sample;
for k=1:K
    shape_est.sample(k,:) = shape_sample(k,:);
    shape_est.mmse(k)     = mean(shape_sample(k,n0+1:end));
    scale_est.sample(k,:) = scale_sample(k,:);
    scale_est.mmse(k)     = mean(scale_sample(k,n0+1:end));   
end

%%% Store results in a structure
if refx>0
    results.psnr = psnr_vec; % peak signal-to-noise ratio
end
if bool==1
    % mean square jump and mean square jump per second
    results.MSJ = [sqrt(MSJ/(Nmc-1)) sqrt(MSJ/(Nmc-1))/((time_vec(end)-time_vec(n0+1))/(Nmc-1))];
end
results.TRF_neglog   = TRF_neglog_vec;                                 % conditional proba of TRF
results.time_vec     = time_vec;                                       % time for each iteration
results.scale.mmse   = scale_est.mmse.*rescale(1).^shape_est.mmse;     % MMSE of scale params 
results.scale.sample = scale_est.sample.*rescale(1).^shape_est.sample; % scale params for each iteration
results.shape        = shape_est;                                      % shape params (MMSE and for each iteration)
results.sig2.mmse    = rescale(rescale(sig2_est.mmse));                % MMSE of noise variance
results.sig2.sample  = rescale(rescale(sig2_est.sample));              % noise variance for each iteration
results.seg          = remove_padding(z_map);                          % estimated segmentation
results.TRF          = rescale(remove_padding(x_hat));                 % estimated TRF
results.n0           = n0;                                             % number of burn-in iterations
    
%%%%%%%%%%%%%%%
function res = compute_TRF_neglog(x,z,sh,sc,var)
    %computes the conditional probability of TRF
    x_small = remove_padding(x);
    z       = remove_padding(z);
    res = 0;
    for l=1:K
      res = res + sum(abs(x_small(z==l)).^sh(l))/sc(l);
    end
    res = res + ...
        0.5*norm(remove_padding(y)-remove_padding(real(ifft2(DFT.*fft2(x)))),'fro')^2/var;
end

end