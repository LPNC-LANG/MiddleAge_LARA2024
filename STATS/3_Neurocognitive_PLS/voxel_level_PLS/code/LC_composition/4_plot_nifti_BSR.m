% Clément Guichet, UGA CNRS UMR 5105 LPNC, 2023

clc
clearvars

% Change accordingly
% myPLS_contrastBehavInteract_norm1-1_LC1_BSR_saliences
% myPLS_contrastBehavInteract_norm1-1_LC2_BSR_saliences

% LC1_BSR_negative_weighted_X % For the clusters of LC1
% LC2_BSR_negative_weighted_X % For the clusters of LC2

% Example
file_name_weighted = 'myPLS_contrastBehavInteract_norm1-1_LC1_BSR_saliences'; 
load('myobj_axial.mat');
slices=-65:3:-30; % Change accordingly

%   Load loadings volume to get min/max values
    clear Ai A myMin myMax
    Ai = spm_vol([file_name_weighted, '.nii']);
    A = spm_read_vols(Ai);
    this_max = max(A(:));
    this_min = min(A(:));
    
    S = spm_read_vols(spm_vol('E:/Research_Projects/MiddleAge_LARA/STATS/3_Neurocognitive_PLS/voxel_level_PLS/code/PLS_analysis/masks/TW_FA_Gaussian25_155subj_mean.nii'));
    s_max=max(S(:));
    
%  Overlay loadings on template brain
    figure;
    set(gcf,'name',file_name_weighted);
    set(gcf,'Position',[484 77 560 751]);
    myobj.img(1).vol = spm_vol('E:/Research_Projects/MiddleAge_LARA/STATS/3_Neurocognitive_PLS/voxel_level_PLS/code/PLS_analysis/masks/TW_FA_Gaussian25_155subj_mean.nii'); % template
    myobj.img(2).vol = Ai; % loadings' positive values
    myobj.img(3).vol = Ai; % loadings' negative values
    myobj.img(1).range=[0 0.5*s_max];
    myobj.img(2).range=[3 5];
    myobj.img(3).range=[-3 -7];
    myobj.figure = gcf;
    myobj.slices=slices;
    myobj.xslices=4;
    paint(myobj);
   