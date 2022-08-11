#!/bin/tcsh
setenv DATA /autofs/cluster/iaslab/FSMAP/FSMAP_data
setenv SCRIPTPATH /autofs/cluster/iaslab/users/danlei/FSMAP/scripts
setenv IMAGE /autofs/cluster/iaslab/users/jtheriault/singularity_images/jtnipyutil/jtnipyutil-2019-01-03-4cecb89cb1d9.simg
setenv PROJNAME wm_nback_block_contrast
setenv SINGULARITY /usr/bin/singularity
setenv SUBJ sub-$1
setenv OUTPUT /autofs/cluster/iaslab/FSMAP2/FSMAP_data

mkdir -p $OUTPUT/BIDS_modeled/$PROJNAME
# mkdir -p /scratch/$USER/$SUBJ/$PROJNAME/BIDS_preproc/$SUBJ/anat
# mkdir -p /scratch/$USER/$SUBJ/$PROJNAME/BIDS_preproc/$SUBJ/func
mkdir -p /scratch/$USER/$SUBJ/$PROJNAME/wrkdir/
# mkdir -p /scratch/$USER/$SUBJ/$PROJNAME/BIDS_modeled

# rsync -ra $DATA/BIDS_fmriprep/fmriprep/ses-01/$SUBJ/func/*wm* /scratch/$USER/$SUBJ/$PROJNAME/BIDS_preproc/$SUBJ/func
# rsync -ra $DATA/BIDS_fmriprep/fmriprep/ses-01/$SUBJ/anat/sub-*_T1w_space-MNI* /scratch/$USER/$SUBJ/$PROJNAME/BIDS_preproc/$SUBJ/anat
# rsync -r $VALIDITYPATH /scratch/$USER/$SUBJ/$PROJNAME/wrkdir/

rsync $SCRIPTPATH/model/$PROJNAME/{wm_lvl1_model.py,wm_lvl1_model_startup.sh} /scratch/$USER/$SUBJ/$PROJNAME/wrkdir/
chmod +x /scratch/$USER/$SUBJ/$PROJNAME/wrkdir/wm_lvl1_model_startup.sh
cd /scratch/$USER
chmod a+rwx /autofs/cluster/iaslab2/FSMAP/FSMAP_data/BIDS_modeled/$PROJNAME

$SINGULARITY exec  \
--bind "$DATA/BIDS_fmriprep/fmriprep/ses-01:/scratch/data" \
--bind "/autofs/cluster/iaslab2/FSMAP/FSMAP_data/xtra_confounds_wm/aCompCor:/scratch/data/aCompCor" \
--bind "$OUTPUT/BIDS_modeled:/scratch/output" \
--bind "/scratch/$USER/$SUBJ/$PROJNAME/wrkdir:/scratch/wrkdir" \
$IMAGE \
/scratch/wrkdir/wm_lvl1_model_startup.sh

# rsync -r /scratch/$USER/$SUBJ/$PROJNAME/BIDS_modeled/ /autofs/cluster/iaslab2/FSMAP/FSMAP_data/BIDS_modeled/
chmod a+rwx /autofs/cluster/iaslab2/FSMAP/FSMAP_data/BIDS_modeled/$PROJNAME

rm -r /scratch/$USER/$SUBJ/$PROJNAME/
exit

# scp -r /scratch/dz609/sub-127/painAvd_CSUS1snegneu/BIDS_preproc/painAvd_CSUS1snegneu/full_model_wf/lvl_one_pipe/modelfit/_subject_id_sub-127_task-pain3_run-05/modelgen/mapflow/_modelgen0 /autofs/cluster/iaslab/users/danlei/temp/
# rsync -r "dz609@door.nmr.mgh.harvard.edu:/autofs/cluster/iaslab/FSMAP/FSMAP_data/BIDS_fmriprep/fmriprep/ses-01/sub-030/func/sub-030_task-wm_events.tsv" "/Users/chendanlei/Desktop/"
