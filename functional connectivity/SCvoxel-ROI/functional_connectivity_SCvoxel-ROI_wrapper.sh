#!/bin/tcsh
setenv DATA /autofs/cluster/iaslab/FSMAP2/FSMAP_data/BIDS_modeled/wm_nback_block/SC_warp/res4d/
setenv SCRIPTPATH /autofs/cluster/iaslab/users/danlei/FSMAP/scripts
setenv SCPATH /autofs/cluster/iaslab/FSMAP2/FSMAP_data/BIDS_modeled/subject_SC_mask/wm
setenv ROIPATH /autofs/cluster/iaslab/users/danlei/FSMAP/scripts/model
setenv IMAGE /autofs/cluster/iaslab/users/jtheriault/singularity_images/jtnipyutil/jtnipyutil-2019-01-03-4cecb89cb1d9.simg
setenv SINGULARITY /usr/bin/singularity
setenv PROJNAME wm_conn_SCvoxel-ROI
setenv OUTPUT /autofs/cluster/iaslab/FSMAP2/FSMAP_data/BIDS_modeled/wm_nback_block/connectivity
setenv ROI $1
# rois = ['FEF', 'LIP', 'LGN', 'V1', 'MGN', 'A1', 'LIPd', 'LIPv']

mkdir -p /scratch/$USER/$PROJNAME/wrkdir/
mkdir -p /scratch/$USER/$PROJNAME/output/
mkdir -p /scratch/$USER/$PROJNAME/data
mkdir -p $OUTPUT

# rsync /autofs/cluster/iaslab/users/danlei/roi/WholeBrain/mni_icbm152_gm_tal_nlin_asym_09a_thresh10.nii /scratch/$USER/$PROJNAME/wrkdir/

# rsync -ra $SCRIPTPATH/model/search_region.nii /scratch/$USER/$PROJNAME/wrkdir/
rsync /autofs/cluster/iaslab/users/danlei/FSMAP/scripts/model/$PROJNAME/* /scratch/$USER/$PROJNAME/wrkdir/
chmod a+rwx /scratch/$USER/$PROJNAME/wrkdir/functional_connectivity_SCvoxel-ROI_startup.sh
cd /scratch/$USER

$SINGULARITY exec \
--bind "$DATA/warped_files:/scratch/data/res4d" \
--bind "$SCPATH/avg_template:/scratch/data/SC_mask" \
--bind "$ROIPATH/wm_conn_roi:/scratch/data/ROI_mask" \
--bind "/autofs/cluster/iaslab/FSMAP/FSMAP_data/BIDS_fmriprep/fmriprep/ses-01/:/scratch/data/confound_files" \
--bind "/scratch/$USER/$PROJNAME/wrkdir:/scratch/wrkdir" \
--bind "$OUTPUT/SCvoxel-ROI2:/scratch/output" \
$IMAGE\
/scratch/wrkdir/functional_connectivity_SCvoxel-ROI_startup.sh

# rsync -r /scratch/$USER/$PROJNAME/output/* /autofs/cluster/iaslab2/FSMAP/FSMAP_data/BIDS_modeled/$PROJNAME/connectivity_analyses/

# cd /autofs/cluster/iaslab/users/danlei/FSMAP/scripts/model/
# chmod -R a+rwx *

rm -r /scratch/$USER/
exit


# scp -r * /autofs/cluster/iaslab/users/danlei/test/
# scp dz609@door.nmr.mgh.harvard.edu:/autofs/cluster/iaslab/users/danlei/test/sub-014_emo3_combined_trial_copes_neg_smoothed3mm_masked_subject_SC_mask.npy /Users/chendanlei/Desktop/
