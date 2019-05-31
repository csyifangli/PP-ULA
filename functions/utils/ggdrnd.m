function x = ggdrnd(nr,nc,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference
% ''A Pratical procedure tu estimate the shape parameter in the generalized
% gaussian distribution''.
% Data are arranged in a [Dim(1) Dim(2)] matrix.
%
% INPUT
% nr = number of rows
% nc = number of columns
% m = mean value (opt, default m = 0)
% scale = scale parameter (opt, default var = 1)
% p = shape parameter (opt, default p = 2)
%
% OUTPUT
% x = samples
%
% EXAMPLE OF USAGE
% x = ggdrnd(1,1024,0,1,0.5);
% built a 1 x 1024 vector filled with samples from a GGD distribution with
% zero mean, scale=1 and shape = 0.5.
%
% Original: Martino Alessandrini, 
% Update: Ningning ZHAO
% Date: Oct 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input arguments
if (nargin < 3) || isempty(varargin{1}), ...
    m = 0; else m = varargin{1}; end
if (nargin < 4) || isempty(varargin{2}), ...
    scale = 1; else scale = varargin{2}; end
if (nargin < 5) || isempty(varargin{3}), ...
    p = 2; else p = varargin{3}; end
%% samples from gamma distribution
% parameters for gamma distribution
lam = 1/p;
z = gamrnd(lam,scale,nr,nc);
%% GGD
w = z.^lam;
y = (-1).^binornd(1,0.5,nr,nc).*w;
x = m + y;
end