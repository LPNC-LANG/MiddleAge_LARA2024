#!/bin/bash

################################################################
# MUST FOLLOW BIDS ARCHITECTURE:
# sub
#   -anat
#       -*T1w.nii.gz
#   -dwi
#       -*.bvec
#       -*.bval
#       -*.json
#       -*dwi.nii.gz
################################################################

NPROC=$(nproc)
############################### STEP 1 ###############################
#             Convert data to .mif format and denoise                #
######################################################################

# Also consider doing Gibbs denoising (using mrdegibbs). Check your diffusion data for ringing artifacts before deciding whether to use it
for_each -nthreads $NPROC -info */dwi : mrconvert IN/*dwi.nii.gz IN/dwi.mif 
for_each -nthreads $NPROC -info */dwi : mrconvert IN/dwi.mif -fslgrad IN/*.bvec IN/*.bval IN/dwi_header.mif 

for_each -nthreads $NPROC -info */dwi : dwidenoise IN/dwi_header.mif IN/dwi_den.mif -noise IN/noise.mif 
for_each -nthreads $NPROC -info */dwi : mrdegibbs IN/dwi_den.mif IN/dwi_den_unr.mif 

# Extract the b0 images from the diffusion data acquired in the AP direction
for_each -nthreads $NPROC -info */dwi : dwiextract IN/dwi_den.mif - -bzero \| mrmath - mean IN/mean_b0_AP.mif -axis 3 

######################################################################
# Runs the dwipreproc command, which is a wrapper for eddy and topup.
#### !!! Here the CAMCAN dataset does not provide reverse encoding, hence -rpe_none !!! ####
#### !!! $NPROC/8 means you should divide by 8 as each subject's preprocessing will be performed by 8 threads already (see at the end of the line) #### !!!
######################################################################
for_each 8 -info */dwi : dwifslpreproc IN/dwi_den.mif IN/dwi_den_preproc.mif -pe_dir AP -rpe_none -readout_time 0.0342002 -eddy_options " --slm=linear --data_is_shelled"  -nthreads 8

# Performs bias field correction. Needs ANTs to be installed in order to use the "ants" option (use "fsl" otherwise)
for_each -nthreads $NPROC -info */dwi : dwibiascorrect ants IN/dwi_den_preproc.mif IN/dwi_den_preproc_unbiased.mif -bias IN/bias.mif 

########################### STEP 2 ###################################
#             Basis function for each tissue type                    #
######################################################################

#The "dhollander" function is best used for multi-shell acquisitions; it will estimate different basis functions for each tissue type. For single-shell acquisition, use the "tournier" function instead
for_each -nthreads $NPROC -info */dwi : dwi2response dhollander IN/dwi_den_preproc_unbiased.mif IN/wm.txt IN/gm.txt IN/csf.txt -voxels IN/voxels.mif 

# Create an average basis function from the subject's DWI data. 
responsemean */dwi/wm.txt ./group_average_wm.txt
responsemean */dwi/gm.txt ./group_average_gm.txt
responsemean */dwi/csf.txt ./group_average_csf.txt

#Upsample the difusion image for better resolution and tracto later
for_each -nthreads $NPROC -info */dwi : mrgrid IN/*unbiased.mif regrid -vox 1.5 IN/dwi_unbiased_upsampled.mif

# Create a mask for future processing steps
for_each -nthreads $NPROC -info */dwi : dwi2mask IN/*unbiased_upsampled.mif IN/mask_up.mif

# Performs multishell-multitissue constrained spherical deconvolution, using the basis functions estimated above
for_each -nthreads $NPROC -info */dwi : dwi2fod msmt_csd IN/*unbiased_upsampled.mif -mask IN/mask_up.mif group_average_wm.txt IN/wmfod_up.mif group_average_gm.txt IN/gmfod_up.mif group_average_csf.txt IN/csffod_up.mif 

# Creates an image of the fiber orientation densities overlaid onto the estimated tissues (Blue=WM; Green=GM; Red=CSF)
# You should see FOD's mostly within the white matter. These can be viewed later with the command "mrview vf.mif -odf.load_sh wmfod.mif"
for_each -nthreads $NPROC -info */dwi : mrconvert -coord 3 0 IN/wmfod_up.mif - \| mrcat IN/csffod_up.mif IN/gmfod_up.mif - IN/vf_up.mif 

# Now normalize the FODs to enable comparison between subjects
for_each -nthreads $NPROC -info */dwi : mtnormalise IN/wmfod_up.mif IN/wmfod_norm_up.mif IN/gmfod_up.mif IN/gmfod_norm_up.mif IN/csffod_up.mif IN/csffod_norm_up.mif -mask IN/mask_up.mif 

########################### STEP 3 ###################################
#            Create a GM/WM boundary for seed analysis               #
######################################################################

# Convert the anatomical image to .mif format, and then extract all five tissue catagories (1=GM; 2=Subcortical GM; 3=WM; 4=CSF; 5=Pathological tissue)
for_each -nthreads 8 -info */dwi : mrconvert IN/../anat/*T1w.nii.gz IN/T1.mif 
for_each -nthreads $NPROC/8 -info */dwi : 5ttgen fsl IN/T1.mif IN/5tt_nocoreg.mif -nthreads 8 
for_each -nthreads $NPROC -info */dwi : mrconvert IN/5tt_nocoreg.mif IN/5tt_nocoreg.nii.gz

# The following series of commands will take the average of the b0 images (which have the best contrast), convert them to NIFTI format, and use it for coregistration.
for_each -nthreads $NPROC -info */dwi : dwiextract IN/*unbiased_upsampled.mif - -bzero \| mrmath - mean IN/mean_b0_processed_up.mif -axis 3 
for_each -nthreads $NPROC -info */dwi : mrconvert IN/mean_b0_processed_up.mif IN/mean_b0_processed_up.nii.gz 


# Uses FSL commands fslroi and flirt to create a transformation matrix for regisitration between the tissue map and the b0 images
for_each -nthreads $NPROC -info */dwi : fslroi IN/5tt_nocoreg.nii.gz IN/5tt_vol0.nii.gz 0 1 #Extract the first volume of the 5tt dataset (since flirt can only use 3D images, not 4D images)
for_each -nthreads $NPROC -info */dwi : flirt -in IN/mean_b0_processed_up.nii.gz -ref IN/5tt_vol0.nii.gz -interp nearestneighbour -dof 6 -omat IN/diff2struct_fsl_up.mat
for_each -nthreads $NPROC -info */dwi : transformconvert IN/diff2struct_fsl_up.mat IN/mean_b0_processed_up.nii.gz IN/5tt_nocoreg.nii.gz flirt_import IN/diff2struct_mrtrix_up.txt 
for_each -nthreads $NPROC -info */dwi : mrtransform IN/5tt_nocoreg.mif -linear IN/diff2struct_mrtrix_up.txt -inverse IN/5tt_coreg_up.mif 

#Create a seed region along the GM/WM boundary
for_each -nthreads $NPROC -info */dwi : 5tt2gmwmi IN/5tt_coreg_up.mif IN/gmwmSeed_coreg_up.mif


########################## STEP 4 ###################################
#                 Run the streamline analysis                        #
######################################################################

# MRtrix3 recommend about 100 million tracks. Here I use 10 million, if only to save time. Read their papers and then make a decision
for_each -nthreads 8 -info */dwi : tckgen -act IN/5tt_coreg_up.mif -backtrack -seed_gmwmi IN/gmwmSeed_coreg_up.mif -nthreads 8 -maxlength 250 -cutoff 0.06 -select 10000k IN/wmfod_norm_up.mif IN/tracks_10M_up.tck 

# Extract a subset of tracks (here, 200 thousand) for ease of visualization
# tckedit tracks_10M.tck -number 200k smallerTracks_200k.tck

# Reduce the number of streamlines with tcksift
for_each -nthreads 8 -info */dwi : tcksift2 -act IN/5tt_coreg_up.mif -out_mu IN/sift_mu_up.txt -out_coeffs IN/sift_coeffs_up.txt -nthreads 8 IN/tracks_10M_up.tck IN/wmfod_norm_up.mif IN/sift_1M_up.txt 