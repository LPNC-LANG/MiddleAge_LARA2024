#!/bin/bash

#################################################################
# 1 Subject-specific WM FOD map's coregistration to template space
#################################################################
# Path to preprocessed data
cd /mnt/h//CAMCAN/DTI/DTI_preproc <-- takes wmfod_norm_up & mask_up as input

for_each sub* : mrregister -transformed /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/2_FOD_template_registration/output/transformed_template/IN_transformed.mif.gz 
	-nl_warp_full /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/2_FOD_template_registration/output/warped_template/1_warp/IN_warp.mif.gz IN/dwi/wmfod_norm_up.mif 
	-mask1 IN/dwi/mask_up.mif /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/1_Generate_FOD_template/output/wmfod_norm_up_template.mif

#################################################################
# 2 Register tracks to template space with FOD-computed warps
#################################################################
for_each * : warpconvert PRE.mif.gz warpfull2deformation -template /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/1_Generate_FOD_template/output/wmfod_norm_up_template.mif
	/mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/2_FOD_template_registration/output/warped_template/2_warpfull2deformation/PRE_warp2fulldeformation.mif.gz
for_each * : warpinvert /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/2_FOD_template_registration/output/2_warpfull2deformation/PRE_warpfull2deformation.mif.gz 
	/mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/2_FOD_template_registration/output/warped_template/3_warped_invert/PRE_warpfull2deformation_invert.mif.gz

# Path to preprocessed data
cd /mnt/h//CAMCAN/DTI/DTI_preproc <-- takes tractograms in subject space as input
for_each * : tcktransform IN/dwi/tracks_10M_up.tck /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/2_FOD_template_registration/output/warped_template/3_warped_invert/PRE_warp_warp2fulldeformation_invert.mif.gz 
	/mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/2_FOD_template_registration/output/Tractograms_in_template_space/PRE_tract_in_template_space.tck

#################################################################
# 3 Generate & Register tensor-based images to template space with FOD-computed warps
#################################################################

# Path to preprocessed data
cd /mnt/h//CAMCAN/DTI/DTI_preproc
for_each -info sub* : dwi2tensor IN/dwi/*unbiased_upsampled.mif - \| tensor2metric - -fa - \| mrcalc - -abs /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/2_FOD_template_registration/data/FA_map_in_subject_space/IN_FA_map.mif.gz


cd /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/2_FOD_template_registration/data/FA_map_in_subject_space
for sub in `cat subjList.txt`; do
 	mrtransform -template /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/1_Generate_FOD_template/output/wmfod_norm_up_template.mif 
 		-warp /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/2_FOD_template_registration/output/warped_template/2_warpfull2deformation/${sub}_warp_warp2fulldeformation.mif.gz 
 		${sub}_FA_map.mif.gz /mnt/e/Research_Projects/MiddleAge_LARA/DWI/TWI/2_FOD_template_registration/output/FA_in_template_space/${sub}_FA_template.nii.gz
done