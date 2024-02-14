#!/bin/bash

#################################################################
# TW-FA images
#################################################################

for sub in `cat subjList.txt`; do
  # SIFT2 weights are the same if generated in subject space or in template space after having applied FOD-computed warps
  tckmap /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/2_FOD_template_registration/output/Tractograms_in_template_space/${sub}*.tck 
  	-info 
  	-tck_weights_in /mnt/h//CAMCAN/DTI/DTI_preproc/${sub}/dwi/sift_1M_up.txt 
  	-stat_tck gaussian -fwhm_tck 25 
  	-template /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/1_Generate_FOD_template/output/wmfod_norm_up_template.mif 
    -contrast scalar_map 
    -image /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/2_FOD_template_registration/output/FA_in_template_space/${sub}_FA_template.nii.gz 
    -stat_vox mean 
    # (optionally)
    # -precise 
    /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/3_TWFA/output/${sub}_TW_FA_Gaussian25.nii.gz
done

# Concatenate volumes along axes 3 (4th dim)
# does not work in input files are in mif.gz
mrcat * TW_FA_Gaussian25_155subj.mif.gz -axis 3