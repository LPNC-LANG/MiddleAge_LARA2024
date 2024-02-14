PLS PIPELINE

#######################
Voxel-level analysis:
#######################

Run (PLSC_main.m)
Save diagnostic plots as svg files
Save res.mat & save_opts.mat
Store Lx and Ly in latent_scores.csv to plot Figures (plot_latent_scores.R)
Report  bootstrap ratios and plot them (plot_BSR_scores.R)


# MRI postproc ----
(salience_mask.m)
Create a salience mask of the negative white matter LC1 alterations and run the cluster analysis 
Identify the number of significantly large clusters (>~1000 voxels)
Generate a mask out of each cluster

(salience_track_overlay.sh) & (salience_track_overlay.R) 
Determine the composition of each cluster
Run (plot_salience_mask.m) to plot figures


#######################
Network-level analysis:
#######################

IDEM

# MRI postproc ----
Reconstruct salience map (salience_recon.m)