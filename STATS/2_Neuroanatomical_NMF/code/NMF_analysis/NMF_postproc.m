% Clément Guichet, UGA CNRS UMR 5105 LPNC, 2023

%% Write one nifti volume for each network
clear all
clc

load('E:\Research_Projects\MiddleAge_LARA\STATS\2_Neuroanatomical_NMF\output\NMF_results\W.mat');
W = W(1,8); % Get the 16-part solution
W = W{1};

% Replace by the row index of the maximum value
for i = 1:16
    [~,row_idx] = max(W');
    row_idx(row_idx~=i) = 0;
    mask_file = 'E:\Research_Projects\MiddleAge_LARA\STATS\3_Neurocognitive_PLS\MRI\TW_FA_Gaussian25_155subj_mean_mask_95.nii';

    % Load mask
    mask_hdr=spm_vol(mask_file);
    mask = spm_read_vols(mask_hdr);
    mask_idx = logical(mask);

    % Binary volume
    file_name = sprintf('covariance_networks_%d', i);
    % Constrain volume to mask
    var_3D = zeros(size(mask));
    var_3D(mask_idx) = row_idx;
    Vi = mask_hdr;
    Vi.dt = [spm_type('float32') 0];
    Vi.fname = [file_name '.nii'];
    spm_write_vol(Vi,var_3D);
    
    % Weighted volume with raw NMF loadings
    file_name_weighted = sprintf('covariance_networks_weighted_%d', i);
    
    var_3D_weighted = zeros(size(mask));
    var_3D_weighted(mask_idx) = W(:,i);
    Vi = mask_hdr;
    Vi.dt = [spm_type('float32') 0];
    Vi.fname = [file_name_weighted '.nii'];
    spm_write_vol(Vi,var_3D_weighted);
end


%% Create a parcellation
clear all
clc

load('E:\Research_Projects\MiddleAge_LARA\STATS\2_Neuroanatomical_NMF\output\NMF_results\W.mat');
W = W(1,8);
fac_matrix = W{1};

[~,row_idx] = max(fac_matrix'); % Discretize probabilities into networks

mask_file = 'E:\Research_Projects\MiddleAge_LARA\STATS\3_Neurocognitive_PLS\MRI\TW_FA_Gaussian25_155subj_mean_mask_95.nii';

% Load mask
mask_hdr=spm_vol(mask_file);
mask = spm_read_vols(mask_hdr);
mask_idx = logical(mask);
       
file_name = ['Network_Parcellation_NMF'];
% Constrain parcellation to mask
var_3D = zeros(size(mask));
var_3D(mask_idx) = row_idx;
Vi = mask_hdr;
Vi.dt = [spm_type('float32') 0];
Vi.fname = [file_name '.nii'];
spm_write_vol(Vi,var_3D);

%% Repeat the two last sections while considering N voxel clusters only within the network
clear all
clc

% Replace by the row index of the maximum value
for network_value = 1:16
    
    vol = spm_vol(sprintf('covariance_networks_%d.nii', network_value));
    [Y, ~] = spm_read_vols(vol); % load a given network
    
    CC = bwconncomp(Y); % find the connected components
    stats = regionprops3(CC,"ConvexVolume","VoxelIdxList","VoxelList","Image", 'SurfaceArea');
    for i = 1:size(stats,1)
        stats{i,6} = size(stats{i,2}{1,1},1);
    end
    stats = sortrows(stats,"Var6", 'descend'); % sort by the size of the cluster
    
%     nb_voxels_to_keep = nnz(Y) * 0.5; % keep clusters until we reach 80% of the voxels of the network
%     sum_voxels = cumsum(nb_current_row);
%     k = find(sum_voxels >= nb_voxels_to_keep, 1);

    for i = 1:(size(stats,1)) % loop through the column containing the voxel list
        % Relative thresholding
%         if stats{i,6} / max(stats{:,6}) >= .20 % keep only the clusters with size which are at least 20% of the maximum one
%             k = i;
%         end
%       Absolute thresholding
        if size(stats{i,2}{1,1},1) >= 25 % keep only the clusters with at minimum N voxels
            k = i;
        end
    end
    
    LCC = [];
    for j = 1:k
        tmp = stats{j,2}{1,1};
        LCC = cat(1, LCC, tmp); % concatenate the voxel indices
    end
    
    clear template_bin
    mask_idx = logical(Y); % binarize
    vectorized_mask = mask_idx(:)';
    size_mask = 1:size(vectorized_mask,2);
    mask_indices = ismember(size_mask, LCC(:,1)); % check if each voxel is part of the input X0 matrix
    template_bin = double(reshape(mask_indices, [121 152 121])); % reshape into 3D volume
    
    template_bin(template_bin==1) = network_value; % label back the network
    
    Vi = vol;
    Vi.dt = [spm_type('float32') 0];
    Vi.fname = sprintf('clusters_network_%d.nii', network_value);
    spm_write_vol(Vi, template_bin); % write with the original hdr input
end

%% Create NMF parcellation with these clusters

% initialize by creating a discrete pseudo-W matrix
pseudo_fac_matrix = zeros(92287,16); % size of the mask

% Load mask
mask_file = 'E:\Research_Projects\MiddleAge_LARA\STATS\3_Neurocognitive_PLS\MRI\TW_FA_Gaussian25_155subj_mean_mask_95.nii';
mask_hdr=spm_vol(mask_file);
mask = spm_read_vols(mask_hdr);
mask_idx = logical(mask);

for network_value = 1:16
    vol = spm_vol(sprintf('clusters_network_%d.nii', network_value));
    [Y, ~] = spm_read_vols(vol); % load a given clustered network
    vec_Y = Y(:); % vectorize for concatenation
    vec_Y_masked = vec_Y(mask_idx(:));
    pseudo_fac_matrix(:,network_value) = vec_Y_masked; 
end

row_idx = sum(pseudo_fac_matrix,2); % Discretize into networks
   
file_name = ['Network_Parcellation_NNMF_cluster'];
% Constrain parcellation to mask
var_3D = zeros(size(mask));
var_3D(mask_idx) = row_idx;
Vi = mask_hdr;
Vi.dt = [spm_type('float32') 0];
Vi.fname = [file_name '.nii'];
spm_write_vol(Vi,var_3D);

