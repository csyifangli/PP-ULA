function [z_new,ind] = label_sample_MRF_4neighbours( z_old,x,K,granu,shape,scale,ch)
%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
% 
% This code has been provided by Ningning Zhao (updated 2014/06/20) and 
% commented by Marie-Caroline Corbineau 05/2019.
%
% Intput: z_old: segmentation labels at previous iteration
%         x: TRF sample
%         K: number of regions
%         granu: granularity parameter for the Potts model
%         shape: GGD shape parameters
%         scale: GGD scale parameters
%         ch: if ch==1 pixels in blackboard are considered
%             if ch==2 pixels in whiteboard are considered
%
% Output: z_new : labels at current ieration (vector)
%         ind : index of z_new in the matrix
%
% This function samples the labels using a Potts Markov field prior. The 
% number of neighbour elements is 4. Using checkeboard, half of the labels 
% can be sampled at a time.
%====================================================================

[m,n]        = size(z_old);
[b{1},b{2}]  = checkerboard(m,n);
pi           = zeros(m,n,K);
board        = zeros(m,n,K);
posterior_pi = zeros(m,n,K);
    
if ch==1
    ind=find(b{1}==1);
elseif ch==2
    ind=find(b{2}==1);
end

for k=1:K
    board(:,:,k)  = neighbours(z_old,k,ch,b{1},b{2});
    if isreal(x)
        pi(:,:,k) = exp(granu*board(:,:,k)).*exp(-abs(x).^shape(k)/scale(k))/(2*scale(k)^(1/shape(k))*gamma(1+1/shape(k)));%GGD (likelihood)
    else
        pi(:,:,k) = exp(granu*board(:,:,k)).*exp(-(abs(real(x)).^shape(k)+abs(imag(x)).^shape(k))/scale(k))/(2*scale(k)^(1/shape(k))*gamma(1+1/shape(k)));%GGD (likelihood)
    end
end
sum_pi = cumsum(pi,3);
sum_pi = sum_pi(:,:,K);
for k=1:K
    posterior_pi(:,:,k) = pi(:,:,k)./sum_pi.*b{ch};
end
posterior_pi = cumsum(posterior_pi,3);
temp         = z_old(ind);
bound        = zeros(numel(temp),K+1);
rho          = rand(numel(temp),1);

bound(:,1) = zeros(size(temp));
temp       = (bound(:,1)<=rho);
for i=2:K+1
    tt         = posterior_pi(:,:,i-1);
    bound(:,i) = tt(ind);
    temp       = temp+(bound(:,i)<=rho);
end
z_new = temp;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b1,b2]=checkerboard(m,n)
N=max(m,n);
pos1          = ones(1,N);
pos1(2:2:end) = 0;
pos1          = toeplitz([pos1]);
pos1=pos1(1:m,1:n);
b1=pos1;
b2=1-b1;
end
%%%
function [board]=neighbours(z,k,ch,b1,b2)

neighbour1 = circshift(z, [0 -1]);
neighbour2 = circshift(z, [0 1]);
neighbour3 = circshift(z, [-1 0]);
neighbour4 = circshift(z, [1 0]);
board=((k-neighbour1)==0) + ((k-neighbour2)==0) + ((k-neighbour3)==0) + ((k-neighbour4)==0);
if ch==1
    board=board.*b1;
elseif ch==2
    board=board.*b2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end





        
    


