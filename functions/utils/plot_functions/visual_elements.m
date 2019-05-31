function [ulc1,ulc2,width,centers,radii,ang] = visual_elements(image)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  image: simulation name
%
% Output: ulc1: upper left hand corner for CNR area 1
%         ulc2: upper left hand corner for CNR area 2
%         width: width of windows for CNR
%         centers: centers of ellipses to display on segmentation
%         radii: radii of ellipses 
%         ang: orientation angle of ellipses
%
% This function outputs coordinates of the two windows used to compute the 
% CNR on the B-mode of the TRF images, and the characteristics of the 
% ellipses highlighting important locations on the segmentation results.
%====================================================================

switch image
    case 'Simu1'
        ulc1     = [120,90];   
        ulc2     = [120,20];   
        width    = [20,20];    
        centers  = [175 155];  
        radii    = [15 50];    
        ang      = [10];       
    case 'Simu2'  
        ulc1     = [85,65];
        ulc2     = [85,145];  
        width    = [22,22];    
        centers  = [200 95; 128 225];
        radii    = [45 20; 20 20];
        ang      = [0;0];
    case 'Kidney'
        ulc1     = [180,310];  
        ulc2     = [180,210]; 
        width    = [20,20];   
        centers  = [113 138];
        radii    = [25 25];
        ang      = [0];
    case 'Thyroid'
        ulc1     = [460,90];    
        ulc2     = [460,128];   
        width    = [40,10];    
        centers = [130 650;115 510];
        radii   = [7 35;6 90];
        ang     = [0;0];
    case 'KidneyReal'
        ulc1     = [213 119];    
        ulc2     = [213 102];   
        width    = [12,7];    
        centers = [83 210];
        radii   = [15 30];
        ang     = [0];
    case 'Bladder'
        ulc1     = [100,15];    
        ulc2     = [100,100];   
        width    = [40,15];    
        centers = [230 135;192 320];
        radii   = [20 60;20 50];
        ang     = [0;60];     
end

