#!/bin/bash

# This code runs with MRtrix3

##########################################################################################################
# 1 Generation of a population template and subject-specific WM FOD map's coregistration to template space
##########################################################################################################

# Create symbolic link to access preprocessed data (fod and mask input data not made available on github)
for sub in `cat template_sample.txt`; do
	echo "${sub}";

	ln -sr /mnt/h/CAMCAN/DTI/DTI_preproc/sub-${sub}/dwi/wmfod_norm_up.mif /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/1_Generate_FOD_template/data/fod_input/${sub}.mif
	ln -sr /mnt/h/CAMCAN/DTI/DTI_preproc/sub-${sub}/dwi/mask_up.mif /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/1_Generate_FOD_template/data/mask_input/${sub}.mif
done 


# Generate population template
population_template /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/1_Generate_FOD_template/data/fod_input -mask_dir /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/1_Generate_FOD_template/data/mask_input /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/1_Generate_FOD_template/output/wmfod_norm_up_template.mif -voxel_size 1.5