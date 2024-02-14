% Clément Guichet, UGA CNRS UMR 5105 LPNC, 2023

clc
clearvars
%% IMPORT Cognitive DATA

% run Import_cog_data.m
% run Import_TIV155.m
AGECOGdata155Subj(:,12) = 1-AGECOGdata155Subj(:,12); % Tip of the tongue ratio: higher is better
AGECOGdata155Subj(:,13) = 1-log(AGECOGdata155Subj(:,13)); % Hotel Task

COG = AGECOGdata155Subj;
%% Import BRAIN DATA - TW_FA Gaussian 25
TW_FA25 = spm_read_vols(spm_vol('E:\Research_Projects\MiddleAge_LARA\DWI\TWI\3_TWFA\output\TW_FA_Gaussian25_155subj.nii'));

for k = 1:size(TW_FA25,4)
    disp(k);
    tmp = TW_FA25(:,:,:,k);
    vec = tmp(:);
    BRAIN_voxels(k,:) = vec;
end

% REDUCE INPUT DIMENSIONALITY of X 
% Filter voxels which are non zeros for most subjects (threshold; 1 = all subjects; .5 = at least half of the subjects)
threshold = .95;
tmp = double(BRAIN_voxels~=0); % binarize
columnMean = mean(tmp,1); % mean for each voxel across subject
BRAIN_voxels_thresholded = BRAIN_voxels(:,columnMean>=threshold); % Input X for PLSC
columnIndices = find(columnMean>=threshold); % Indices in the original matrix indicating the voxel position

%% Create mean binary mask from mean TW-FA maps for plotting
% Files are in the ./masks directory

mask_hdr = spm_vol('TW_FA_Gaussian25_155subj_mean.nii'); % Load mean file
mask = spm_read_vols(mask_hdr);
mask_idx = logical(mask); % binarize
    
vectorized_mask = mask_idx(:)';
size_mask = 1:size(vectorized_mask,2);
mask_indices = ismember(size_mask, columnIndices(1,:)); % check if each voxel is part of the input X0 matrix
template_bin = reshape(mask_indices, [121 152 121]); % reshape into 3D volume

Vi = mask_hdr;
Vi.dt = [spm_type('float32') 0];
Vi.fname = ['TW_FA_Gaussian25_155subj_mean_mask_95.nii'];
spm_write_vol(Vi, template_bin); % write with the original hdr input


%% Create Age groups for contrasts
COG = AGECOGdata155Subj;

for i=1:size(COG,1) % age at time of MRI acquisition
    if COG(i,5)>= 56
        age_groups=3;
    elseif COG(i,5)<=55 && COG(i,5)>=51
        age_groups=2;
    elseif COG(i,5)<51
        age_groups=1;
    end
    grouping(:,i) = age_groups;
end

input.group_names={'45-50',...
                   '51-55',...
                   '56-60'}; 

%% Define all the inputs
% Here we regress out the confounds variables on both modalities, as done
% in Kebets et al. (2019)

confounds = cat(2, COG(:,6:8), participantdataDTIS1); %  handedness gender + -0.5=male/0.5=female + MMSE + TIV

X = BRAIN_voxels_thresholded; % TW_FA Gaussian 25
[X_reg, ~, ~, ~] = CBIG_glm_regress_matrix(X,confounds,0,[]);

Y_age = cat(2, COG(:,[5, 9:16]));
Y_noage = Y_age(:,2:end); % Without age
[reg_Y_noage, ~, ~, ~] = CBIG_glm_regress_matrix(Y_noage,confounds,0,[]);
% quantile normalization to reduce bias due to non-gaussianity of our behavioral variables
qreg_Y_noage = quantilenorm(reg_Y_noage);
% Y_reg = qreg_Y_noage;
Y_reg = cat(2, Y_age(:,1), qreg_Y_noage); % add back age as meta-info

% Grab the labels for visualization
labels_age = {'Age', 'Cattell', 'Proverb', 'Naming', 'ToT_inverse', 'Hotel_Task', 'Sentence', 'Story_rec', 'VF'};
Y_labels = labels_age;
%% Check all inputs for validity
% !!! always run this function to check your setup before running PLS !!!
binary_mask = './masks/TW_FA_Gaussian25_155subj_mean_mask_95.nii';
structural_template = './masks/TW_FA_Gaussian25_155subj_mean.nii';

output_path = '../../output/raw_results'; 

myPLS_inputs_customed_interact
[input,pls_opts,save_opts] = myPLS_initialize(input,pls_opts,save_opts);


%% Run PLS analysis (including permutation testing and bootstrapping)

res = myPLS_analysis(input,pls_opts); 

bootstrap_ratios_cog = res.U./res.boot_results.Ub_std;
res.V = res.V * (-1)/res.boot_results.Vb_std;
myPLS_plot_results(res,save_opts);