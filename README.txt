PP-ULA_joint_deconv_seg_ultrasound

***************************************************************
* author: Marie-Caroline Corbineau    		              *
* institution: CVN, CentraleSupelec, Universite Paris-Saclay  *
* date: Thursday, May 23th 2019   	                      *
* License: GNU General Public License v3.0                    *
***************************************************************

*****************************************************
* RECOMMENDATIONS:                                  *
* This code is designed to work with                *
* Matlab R2018b (or earlier versions)               *
*****************************************************

-----------------------------------------------------------------------------------------------------
DESCRIPTION:

This code allows to perform joint deconvolution and segmentation of radio-frequency ultrasound images. 
In the proposed MCMC approach, the following hybrid Gibbs sampler estimates the segmentation, the noi-
se variance, the deconvolved tissue reflectivity function (TRF), and the shape and scale parameters of
the underlying GGD distributions in the TRF. 

1) Sample the noise variance according to an inverse Gamma distribution;
2) Sample the shape parameters using a Metropolis Hastings random walk;
3) Sample the scale parameter according to an inverse Gamma distribution;
4) Sample the hidden label field following a Potts model;
5) Sample the TRF using PP-ULA.

PP-ULA is an original proximal unadjusted Langevin sampling method. It makes use of a forward-backward
step and of a preconditioner, also called variable metric, whose purpose is to improve the computatio-
nal time. It can be noted that the preconditioner can be set equal to the identity matrix; this confi-
guration is referred to as P-ULA.

This toolbox consists of the following subfolders:
1) data: contains masks and PSF for the simulated images, as well as in vivo ultrasound data
2) functions: contains the sampler programs 
        2.1) evaluation_metrics: contains the functions evaluating the quality of the results
        2.2) Gibbs_sampler_functions: contains the functions for the five steps of the aforementioned 
                                      Gibbs sampler
        2.3) prox_functions: contains functions that are specific to the proposed PP-ULA algorithm
        2.4) utils: contains useful minor functions
             2.4.1) plot_functions: contains functions that are used to present the results
All functions are documented.

Information about the data:
Real in vivo images (Thyroid, Bladder and KidneyReal) are provided as a courtesy by Denis Kouame, 
from Institut de Recherche en Informatique de Toulouse (IRIT).

-----------------------------------------------------------------------------------------------------
SPECIFICATIONS for using PP-ULA_joint_deconv_seg_ultrasound:

A demo file is provided :
* main_script.m runs the proposed method with several configurations (more information within). 
* util_main.m is a program used by main_script.m, this file containes the specifications regarding the 
number of iterations.

-----------------------------------------------------------------------------------------------------
RELATED PUBLICATION:

# M.-C. Corbineau, D. Kouame, E. Chouzenoux, J.-Y. Tourneret, & J.-C. Pesquet (2019). Preconditioned 
P-ULA for Joint Deconvolution-Segmentation of Ultrasound Images. arXiv preprint arXiv:1903.08111.
-----------------------------------------------------------------------------------------------------
