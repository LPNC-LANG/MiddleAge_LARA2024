% Clément Guichet, UGA CNRS UMR 5105 LPNC, Nov 2023

%% IMPORT DATA
clearvars
clc

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

%% Perform Orthogonal Projective NNMF
clc

clear W
clear H

i = 0;
for k = 2:2:20 % granularity
    i = i+1;
    disp(k)
    [W{i}, H{i}] = opnmf_mem(BRAIN_voxels_thresholded', k, NNDSVD(BRAIN_voxels_thresholded',k,3));
end

%% Frobenius norm
close all

load('W.mat')
load('H.mat')

clear sol
for sol = 1:10
    A = BRAIN_voxels_thresholded';
    R = W{sol}*H{sol};
    error(1,sol) = sol*2;
    error(2,sol) =  norm(A - R, 'fro');
    
%     RMSE(1,sol) = sol*2;
%     squaredError = (A-R) .^ 2;
%     meanSquaredError = sum(squaredError(:)) / numel(A);
%     RMSE(2,sol) = sqrt(meanSquaredError);
end

for i = error(1,2:end)
    disp(i) % number of components
    error(3,i/2) = error(2,(i/2)) -  error(2,i/2-1);
end

% Curvature
x = error(1,:);
yy = error(3,:);

figure(1);
fig1_comps.fig = gcf;
hold on
fig1_comps.p1 = plot(x, yy, 'o');
fig1_comps.p2 = plot(x, yy);
hold off
%========================================================
% ADD LABELS, TITLE, LEGEND
title('Recon error of opNMF solutions');
xlabel('K components');
ylabel('F norm');
legendX = .82; legendY = .87; legendWidth = 0.02; legendHeight = 0.02;
fig1_comps.legendPosition = [legendX, legendY, legendWidth, legendHeight];

%% Write W and H
writematrix(H{8}, 'expansion_coeffs.csv');
% fac_matrix = W{8}';

