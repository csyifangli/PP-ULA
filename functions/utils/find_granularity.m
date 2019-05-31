function [granularity] = find_granularity(noise_seed,image,...
    n0,Nmc,precond_bool,threshold,granularity_max,granularity_min)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  noise_seed: noise seed
%         image: simulation name
%         n0: number of burn-in iterations
%         Nmc: number of iterations after burn-in
%         precond_bool: 0 if P-ULA, 1 if PP-ULA
%         threshold: maximal number of isolated point in the segmentation
%         granularity_max: maximal value for the granularity
%         granularity_min: minimal value for the granularity
%
% Output: granularity: final granularity value
%
% This function finds the granularity value that leads to a segmentation
% with a number of isolated points that is close to the threshold.
%====================================================================

max_it        = 20;
points_before = 0;
fprintf('=====> Run PP-ULA several times to find the best granularity parameter for %s\n',image)
for i=1:max_it
    %set the granularity between the maximal and minimal value
    granularity = (granularity_max+granularity_min)/2;
    %run the Hybrid sampler using PP-ULA
    res         = demo(image,n0,Nmc,precond_bool,noise_seed,granularity);
    %compute the number of isolated points in the segmentation
    points      = compute_isolated_points(res.seg);
    fprintf('granularity = %.3e isolated points = %.3e threshold = %.3e\n',granularity,points,threshold)
    if points<threshold
        granularity_max = granularity;
    else
        granularity_min = granularity;
    end
    if i>1 && abs(points_before-points)<1e-4,break,end
    points_before = points;
end

end

