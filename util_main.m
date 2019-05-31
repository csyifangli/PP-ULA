function [] = util_main(mode,path_results,image,varargin)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  mode: value in {1,2,3}
%         path_results: path to the results folder
%         image: simulation name
%         varargin: varargin{1}: method, 'P-ULA' or 'PP-ULA'
%                   varargin{2}: granularity parameter for Potts model
%
% mode = 1, this function plots results.
% mode = 2, this function runs the simulation for the chosen image and method.
% mode = 3, this function runs several times the chosen simulation for
% PP-ULA in order to find the granularity value that leads to a number of
% isolated points in the segmentation image which is close to the given
% threshold.
%====================================================================
addpath(genpath(pwd))
noise_seed = 3;
switch mode
    case 1
        % plot results
        plot_results(image,path_results)
        
    case 2
        method      = varargin{1};
        granularity = varargin{2};
        switch method
            case 'PPULA'
                precond_bool = 1; % with preconditioner 
                switch image
                    case 'Simu1'
                        n0           = 4001; % number of burn-in iterations
                                             % index starts at 1
                        Nmc          = 4000; % number of iterations after burn-in
                    case 'Simu2'
                        n0           = 10001; 
                        Nmc          = 10000;
                    case 'Kidney'
                        n0           = 7001;
                        Nmc          = 7000;
                    case 'Thyroid'
                        n0           = 3001;
                        Nmc          = 3000;
                    case 'Bladder'
                        n0           = 5001;
                        Nmc          = 5000;
                    case 'KidneyReal'
                        n0           = 5001;
                        Nmc          = 5000; 
                end
            case 'PULA'
                switch image
                    case {'Simu1','Simu2'}
                        precond_bool = 0;     % without preconditioner
                        n0           = 70001; 
                        Nmc          = 70000;
                end
        end
        fprintf('=====> Run %s on %s for a total of %d iterations with granularity = %.3f\n',method,image,n0+Nmc-1,granularity)
        % run simulation
        results = demo(image, n0, Nmc, precond_bool,noise_seed,granularity);
        % save results
        save(strcat(path_results,image,'/',...
                'Precond',num2str(precond_bool),...
                '_',image,...
                '_n0_',num2str(n0),'_Nmc_',num2str(Nmc),...
                '_granu_',strrep(num2str(granularity),'.','_')),'results')
 
    case 3 
        precond_bool = 1; % preconditioner
        switch image
            case 'Simu1' 
                n0               = 4001; % number of burn-in iterations
                Nmc              = 4000; % number of iterations after burn-in
                threshold        = 5e-4; % maximal number of isolated points in the segmented image
                granularity_max  = 1.7;  % maximal granularity value
                granularity_min  = 0.2;  % minimal granularity value
            case 'Simu2'
                n0               = 10001; 
                Nmc              = 10000; 
                threshold        = 1e-3;
                granularity_max  = 1.7;
                granularity_min  = 0.2;
            case 'Kidney'
                n0               = 7001;
                Nmc              = 7000; 
                threshold        = 8e-3;
                granularity_max  = 1.7;
                granularity_min  = 0.2;
            case 'Thyroid'
                n0               = 3001; 
                Nmc              = 3000; 
                threshold        = 8e-4; 
                granularity_max  = 2.5;
                granularity_min  = 1.9;
            case 'Bladder'
                n0               = 5001; 
                Nmc              = 5000; 
                threshold        = 8e-4;
                granularity_max  = 2.5;
                granularity_min  = 1.8;
            case 'KidneyReal'
                n0               = 5001; 
                Nmc              = 5000; 
                threshold        = 8e-4;
                granularity_max  = 2.5;
                granularity_min  = 1.8;
        end
        % find and save the granularity value that leads to a number of 
        % isolated points which is close to the threshold
        granularity = find_granularity(noise_seed,image,n0,Nmc,...
            precond_bool,threshold,granularity_max,granularity_min);
        % save the final granularity value
        save(strcat(path_results,image,'/','granularity_',image),'granularity')
    
end
end

