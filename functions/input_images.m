
function [y,rescale,remove_padding,refx,K,DFT,psf_ref,struct] = input_images(image)
%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  image: string, simulation name
%
% Output: y: observed noisy RF image
%         rescale: function that multiplies its input by the scaling factor       
%         remove_padding: function that removes padding from its image input
%         refx: 0: ground-truth is not known, real data
%               1: TRF and segmentation ground-truths are known
%               2: only TRF ground-truth is known
%         K: number of segmentation labels
%         DFT: direct Fourier transform of the PSF
%         psf_ref: points spread function (PSF)
%         struct: structure with various entries depending on refx
%                 struct.cst : constrast to compute B-mode
%                 struct.TRF : true TRF  (if refx>0)        
%                 struct.RF_wo_noise : true RF image without noise (if refx>0)
%                 struct.sig2 : noise variance (if refx=1)
%                 struct.shape : true shape parameters (if refx=1)
%                 struct.scale : true scale parameters (if refx=1)
%                 struct.mask : true segmentation (if refx=1)
%
% This function generates or loads the data corresponding to a simulation.
%====================================================================

image_var = 10; % set the TRF image variance equal to 10
                % to avoid scaling effects

switch image
    case 'Simu1'
        cst     = 0.2; 
        padding = [0,0];
        refx    = 1;
        shape   = [0.6 1.5];
        scale   = [1 1];
        K       = 2;
        load('mask_2regions','mask')
        [m,n] = size(mask);
        TRF   = zeros(m,n);
        % generate generalized Gaussian distribution for every region
        for i=1:K
            temp         = ggdrnd(m,n,0,scale(i),shape(i));
            TRF(mask==i) = temp(mask==i);
        end
    
    case 'Simu2'
        cst     = 5;
        padding = [6,6];
        refx    = 1; 
        shape   = [0.5 1 1.5];
        scale   = [4 50 100];
        K       = 3;
        load('mask_3regions','mask')
        [m,n] = size(mask);
        TRF   = zeros(m,n);
        % generate generalized Gaussian distribution for every region
        for i=1:K
            temp          = ggdrnd(m,n,0,scale(i),shape(i));
            TRF(mask==i)  = temp(mask==i);
        end
        
    case 'Kidney'
        cst     = 1; 
        padding = [6,6];
        refx    = 2; 
        K       = 3; 
        %%% create RF image using code from Field II and optical photo of
        %%% human kidney and liver
        TRF = human_kidney_phantom_TRF(1e6); % 1e6 scattering points
        TRF = TRF(156:449,136:489);
        TRF = TRF-mean2(TRF);
        TRF = TRF./sqrt(mean2(TRF(:).^2))*image_var;
    
    case 'Thyroid'
        cst     = 400; 
        padding = [10,5];
        refx    = 0;
        K       = 3;
        load('Thyroid','image_RF','apsf'); 
        % crop RF image
        y = image_RF(81:950,99:238);
        % estimate psf from RF image of wires
        psf_ref = apsf(745:753,120:134)-mean2(apsf);
            
    case 'Bladder' 
        cst     = 5;
        padding = [40,10];
        refx    = 0;
        K       = 3;
        load('Bladder','image_RF','apsf');
        % crop RF image
        y = image_RF(290:659,:); 
        % estimate psf from RF image of wires
        psf_ref = apsf(620:626,159:161)-mean2(apsf);
        
    case 'KidneyReal'
        cst     = 5; % constrast to plot B-mode
        padding = [20,10];
        refx    = 0;
        K       = 2;
        load('KidneyReal','image_RF','apsf')
        % crop RF image
        y = image_RF(300:649,25:224);
        % estimate psf from RF image of a wire
        psf_ref = apsf(620:626,159:161)-mean2(apsf);
end

if refx>0
    load('psf_simu','psf_ref')
    TRF_temp = padarray(TRF,[padding(1),padding(2)],'symmetric');
    [m,n]    = size(TRF_temp);
else
    y = y-mean2(y);
    if strcmp(image,'Thyroid')
        y = padarray(y,[padding(1),0],0);
        y = padarray(y,[0,padding(2)],'symmetric');
    else
        y = padarray(y,[padding(1),padding(2)],0);
    end
    [m,n] = size(y);
end

% create direct Fourier trasnform of the PSF
% used to compute Hx using fft2
[a,b]   = size(psf_ref);   
DFT     = zeros(m,n);
DFT((m-a+3)/2:(m+a+1)/2,(n-b+3)/2:(n+b+1)/2) = psf_ref;
DFT = fft2(fftshift(DFT)/sum(DFT(:)));

if refx>0
    RF        = real(ifft2(DFT.*fft2(TRF_temp)));
    P_RF      = mean((RF(:)-mean2(RF)).^2);
    BSNR      = 30;
    sig2      = P_RF/10^(BSNR/10);
    % add noise
    y         = RF + sqrt(sig2)*randn(size(RF));
end

% "normalization": introduce a scaling factor to fix the standard
% deviation and avoid being sensitive to scaling effects
% -> estimate x/scale_fact instead of x
scale_fact = sqrt(mean2(y.^2))/image_var;
y          = y./scale_fact;

rescale         = @(x) scale_fact.*x;
remove_padding  = @(x) x(padding(1)+1:end-padding(1),padding(2)+1:end-padding(2));

struct.cst = cst;
if refx>0
    struct.TRF = TRF;        
    struct.RF_wo_noise = RF; 
    if refx==1
        struct.sig2   = sig2;
        struct.shape  = shape;
        struct.scale  = scale;
        struct.mask   = mask;

    end
end

end