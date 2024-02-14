% Clément Guichet, UGA CNRS UMR 5105 LPNC, Oct 2023
clc
clearvars
%% IMPORT DATA

% run Import_cog_data.m
% run Import_TIV155.m
AGECOGdata155Subj(:,12) = 1-AGECOGdata155Subj(:,12); % Tip of the tongue ratio: higher is better
AGECOGdata155Subj(:,13) = 1-log(AGECOGdata155Subj(:,13)); % Hotel Task

COG = AGECOGdata155Subj;

%% Load the subject-specific coefficients for each networks
X = load('E:\Research_Projects\MiddleAge_LARA\STATS\2_Neuroanatomical_NMF\output\NMF_results\H.mat');
X = X.H(1,8);
X = X{1}';
%% Create Age groups for meaningful contrasts
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
output_path = '.'; % change accordingly

myPLS_inputs_customed
[input,pls_opts,save_opts] = myPLS_initialize(input,pls_opts,save_opts);

%% Run PLS analysis (including permutation testing and bootstrapping)

res = myPLS_analysis(input,pls_opts);

bootstrap_ratios_U = res.U./res.boot_results.Ub_std;
bootstrap_ratios_V = res.V./res.boot_results.Vb_std;

