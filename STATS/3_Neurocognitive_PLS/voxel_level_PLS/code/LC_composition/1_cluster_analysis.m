% Clément Guichet, UGA CNRS UMR 5105 LPNC, 2023

%% Create BSR negative masks from LC1
clc 
clearvars

%% Perform cluster analysis
vol = spm_vol('myPLS_contrastBehavInteract_norm1-1_LC1_BSR_saliences.nii');
[Y, ~] = spm_read_vols(vol);

% Modify to put the 56-60 group in the positive values of LC1
Y = Y*(-1);
% Find the salient voxels
Y(Y>-2.58) = 0; % keep negative BSR

CC = bwconncomp(Y,6); % find the connected components
stats = regionprops3(CC,"ConvexVolume","VoxelIdxList","VoxelList","Image", 'SurfaceArea');

% You may need to re-run this block code manually
for i = 1:size(stats,1)
   stats{i,6} = size(stats{i,2}{1,1},1);
   stats{i,7} = stats{i,6}/sum(stats{:,6}); % percentage of total (negative) alterations
end
stats = sortrows(stats,"Var7", 'descend'); % sort by the size of the cluster

%% Generate a binary mask for computing the overlap
for i = 1:3
    LCC = stats{i,2}{1,1};

    clear template_bin
    mask_idx = logical(Y); % binarize
    vectorized_mask = mask_idx(:)';
    size_mask = 1:size(vectorized_mask,2);
    mask_indices = ismember(size_mask, LCC(:,1)); % check if each voxel is part of the input cluster 
    template_bin = double(reshape(mask_indices, [121 152 121])); % reshape into 3D volume

    Vi = vol;
    Vi.dt = [spm_type('float32') 0];
    Vi.fname = sprintf('LC1_BSR_negative_binary_cluster%d.nii', i);
    spm_write_vol(Vi, template_bin); % write with the original hdr input
end

%% Generate a salience map for each cluster
for i = 1:3
    LCC = stats{i,2}{1,1}; % Grab the voxel indices of the cluster
    Y_vec = Y(:); % vectorize
    Y_new = Y_vec(LCC(:,1),:); % keep only the salience values at the indices
    
    % Load mask of the cluster
    mask_file = sprintf('LC1_BSR_negative_binary_cluster%d.nii', i);
    mask_hdr=spm_vol(mask_file);
    mask = spm_read_vols(mask_hdr);
    mask_idx = logical(mask);
    
    var_3D_weighted = zeros(size(mask));
    var_3D_weighted(mask_idx) = Y_new(:,1); % Assign salience values to this mask
    Vi = mask_hdr;
    Vi.dt = [spm_type('float32') 0];
    Vi.fname = sprintf('LC1_BSR_negative_weighted_cluster%d.nii', i);
    spm_write_vol(Vi,var_3D_weighted);
end


%% Create BSR positive masks from LC2
clc 
clearvars

%% Perform cluster analysis
vol = spm_vol('myPLS_contrastBehavInteract_norm1-1_LC2_BSR_saliences.nii');
[Y, ~] = spm_read_vols(vol);

% Find the salient voxels
Y(Y<2.58) = 0; % keep positive BSR

CC = bwconncomp(Y, 6); % find the connected components
stats = regionprops3(CC,"ConvexVolume","VoxelIdxList","VoxelList","Image", 'SurfaceArea');

% You may need to re-run this block code manually
for i = 1:size(stats,1)
   stats{i,6} = size(stats{i,2}{1,1},1);
   stats{i,7} = stats{i,6}/sum(stats{:,6}); % percentage of total (negative) alterations
end
stats = sortrows(stats,"Var7", 'descend'); % sort by the size of the cluster

%% Generate a binary mask for computing the overlap
for i = 1:2
    LCC = stats{i,2}{1,1};

    clear template_bin
    mask_idx = logical(Y); % binarize
    vectorized_mask = mask_idx(:)';
    size_mask = 1:size(vectorized_mask,2);
    mask_indices = ismember(size_mask, LCC(:,1)); % check if each voxel is part of the input cluster 
    template_bin = double(reshape(mask_indices, [121 152 121])); % reshape into 3D volume

    Vi = vol;
    Vi.dt = [spm_type('float32') 0];
    Vi.fname = sprintf('LC2_BSR_positive_binary_cluster%d.nii', i);
    spm_write_vol(Vi, template_bin); % write with the original hdr input
end

%% Generate a salience map for each cluster
for i = 1:2
    LCC = stats{i,2}{1,1}; % Grab the voxel indices of the cluster
    Y_vec = Y(:); % vectorize
    Y_new = Y_vec(LCC(:,1),:); % keep only the salience values at the indices
    
    % Load mask of the cluster
    mask_file = sprintf('LC2_BSR_positive_binary_cluster%d.nii', i);
    mask_hdr=spm_vol(mask_file);
    mask = spm_read_vols(mask_hdr);
    mask_idx = logical(mask);
    
    var_3D_weighted = zeros(size(mask));
    var_3D_weighted(mask_idx) = Y_new(:,1); % Assign salience values to this mask
    Vi = mask_hdr;
    Vi.dt = [spm_type('float32') 0];
    Vi.fname = sprintf('LC2_BSR_positive_weighted_cluster%d.nii', i);
    spm_write_vol(Vi,var_3D_weighted);
end
