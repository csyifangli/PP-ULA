
function [TRF] = human_kidney_phantom_TRF (N)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  N: number of scattering point
%
% Output: TRF: simulated tissue reflectivity function
%
% This function produces a tissue-mimicking phantom with N scatterers that 
% are uniformly distributed over a digital image of human kidney tissue 
% provided with the Field II ultrasound simulator.
% The amplitude of each scatterer is produced using a zero-mean Gaussian 
% distribution whose variance is linked to the amplitude of the point on 
% the digital image.
%
% This code is provided with Field II https://field-ii.dk/
% Slightly modified by M.-C. Corbineau to keep relevant lines 14/12/2018
%====================================================================

% Load the bitmap image
[liv_kid, ~] = bmpread('kidney_cut.bmp');

% Define image coordinates
liv_kid  = liv_kid';
[Nl, Ml] = size(liv_kid);

x_size  = 100/1000 ;    %  Size in x-direction [m]
dx      = x_size/Nl;          %  Sampling interval in x direction [m]
z_size  = 100/1000 ;    %  Size in z-direction [m]
dz      = z_size/Ml;          %  Sampling interval in z direction [m]
y_size  = 15/1000;      %  Size in y-direction [m]
theta   = 35/180*pi;     %  Rotation of the scatterers [rad]
theta   = 0;
z_start = 2/1000;

% Calculate position data
x0 = rand(N,1);
x  = (x0-0.5)* x_size;
z0 = rand(N,1);
z  = z0*z_size+z_start;
y0 = rand(N,1);
y  = (y0-0.5)* y_size; 

%  Find the index for the amplitude value
xindex = round((x + 0.5*x_size)/dx + 1);
zindex = round((z - z_start)/dz + 1);
inside = (0 < xindex)  & (xindex <= Nl) & (0 < zindex)  & (zindex <= Ml);
index  = (xindex + (zindex-1)*Nl).*inside + 1*(1-inside);

% Amplitudes with different variance must be generated according to the the 
% input map.
% The amplitude of the liver-kidney image is used to scale the variance
amp = exp(liv_kid(index)/100);
amp = amp-min(min(amp));
amp = 1e6*amp/max(max(amp));
amp = amp.*randn(N,1).*inside;

TRF        = zeros(size(liv_kid));
TRF(index) = amp;
TRF        = TRF';




