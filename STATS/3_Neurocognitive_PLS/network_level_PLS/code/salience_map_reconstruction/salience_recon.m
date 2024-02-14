%% Create a parcellation from the salient NMF networks (BSR <= -3)
clear all
clc

% Load BSR values for all NMF networks
load('E:\Research_Projects\MiddleAge_LARA\STATS\3_Neurocognitive_PLS\network_level_PLS\output\raw_results\res_NMF.mat');
bootstrap_ratios = res.V./res.boot_results.Vb_std;

% Load NMF solution
load('E:\Research_Projects\MiddleAge_LARA\STATS\2_Neuroanatomical_NMF\output\NMF_results\W.mat');
W = W(1,8);
fac_matrix = W{1};

% Discretize probabilities into networks
[~,row_idx] = max(fac_matrix'); 

% Find the salient NMF networks
bootstrap_ratios([1,2,3,6,8,9,10,12,13,15],:)=0;

% Assign each voxel the corresponding BSR value of its network
for voxel = 1:size(row_idx,2)
    network = row_idx(1,voxel); % grab the network value
    BSR_label(1,voxel) =  bootstrap_ratios(network,1)*-1; % grab the BSR value and assign it to the voxel
end

% Write out the reconstructed salience map
mask_file = 'E:/Research_Projects/MiddleAge_LARA/STATS/3_Neurocognitive_PLS/voxel_level_PLS/code/PLS_analysis/masks/TW_FA_Gaussian25_155subj_mean_mask_95.nii';

% Load mask
mask_hdr=spm_vol(mask_file);
mask = spm_read_vols(mask_hdr);
mask_idx = logical(mask);
       
file_name = ['NNMF_BSR_LC1'];
% Constrain parcellation to mask
var_3D = zeros(size(mask));
var_3D(mask_idx) = BSR_label;
Vi = mask_hdr;
Vi.dt = [spm_type('float32') 0];
Vi.fname = [file_name '.nii'];
spm_write_vol(Vi,var_3D);

%% Plot the reconstruction salience map
load('myobj_axial.mat');
slices=-60:5:15;
    
%   Load loadings volume to get min/max values
    clear Ai A myMin myMax
    Ai = spm_vol([file_name '.nii']);
    A = spm_read_vols(Ai);
    this_max = max(A(:));
    this_min = min(A(:));
    
    S = spm_read_vols(spm_vol('E:/Research_Projects/MiddleAge_LARA/STATS/3_Neurocognitive_PLS/voxel_level_PLS/code/PLS_analysis/masks/TW_FA_Gaussian25_155subj_mean.nii'));
    s_max=max(S(:));

    figure;
    set(gcf,'name',file_name);
    set(gcf,'Position',[484 77 560 751]);
    myobj.img(1).vol = spm_vol('E:/Research_Projects/MiddleAge_LARA/STATS/3_Neurocognitive_PLS/voxel_level_PLS/code/PLS_analysis/masks/TW_FA_Gaussian25_155subj_mean.nii'); % template
    myobj.img(2).vol = spm_vol([file_name '.nii']); % loadings' positive values
    myobj.img(3).vol = spm_vol([file_name '.nii']); % loadings' negative values
    myobj.img(1).range=[0 0.5*s_max];
    myobj.img(2).range=[2.58 3];
    myobj.img(3).range=[-2.58 -7];
    myobj.figure = gcf;
    myobj.slices=slices;
    myobj.xslices=4;
    paint(myobj);