
close all; clear all; clc

%%%%% Modes
% This demo is configured for several modes:
% mode = 1: plots results that are already saved
%        2: runs a method for a fixed granularity parameter
%        3: runs PP-ULA several times to determine the best granularity 
%           parameter given a threshold corresponding to the maximal number
%           of isolated points in the segmentation

%%%%% Experiments
% image = 'Simu1': simulated image, 2 regions
%         'Simu2': simulated image, 3 regions
%         'Kidney': tissue mimicking phantom from Field II
%         'Thyroid': real in vivo data of thyroid flux
%         'Bladder': real in vivo data of a mouse bladder
%         'KidneyReal': real in vivo data of a kidney

%%%%% Methods
% method = 'PPULA': preconditioned proximal unadjusted Langevin algorithm
%          'PULA': proximal unadjusted Langevin algo

% The number of iterations during and after burn-in are specified in
% util_main.m

%% Plot results
mode         = 1;
path_results = 'results/'; % folder where the results are saved
image        = 'Simu1';
if ~exist(strcat(path_results,image), 'dir')
    fprintf('No results are saved for %s\n',image)
else
    util_main(mode,path_results,image)
end

%% Run P-ULA or PP-ULA fr a fix value of the granularity

mode         = 2;
path_results = 'results/'; % folder where the results will be saved
image        = 'Simu1';
method       = 'PPULA';
if ~exist(strcat(path_results,image), 'dir')
    mkdir(strcat(path_results,image))
end

%%% load saved value for granularity
%%% load(strcat('results/',image,'/','granularity_',image))

%%% or fix it by hand
% Below are the granularity values that were used in our paper. 
% These values have been found using mode=3.
 granularity = 1.1141;  % Simu1      (~10min)
% granularity = 1.1023; % Simu2      (~40min)
% granularity = 1.2763; % Kidney     (~40min)
% granularity = 1.9938; % Thyroid    (~35min)
% granularity = 2.4125; % Bladder    (~30min)
% granularity = 2.3250; % KidneyReal (~20min)

util_main(mode,path_results,image,method,granularity)

%% Run PP-ULA several times in order to find the best granularity parameter 
% for the Potts model given a maximal number of isolated points in the
% segmentated image

mode         = 3;
path_results = 'results/';
image        = 'Simu1';
if ~exist(strcat(path_results,image),'dir')
    mkdir(strcat(path_results,image))
end
util_main(mode,path_results,image)
