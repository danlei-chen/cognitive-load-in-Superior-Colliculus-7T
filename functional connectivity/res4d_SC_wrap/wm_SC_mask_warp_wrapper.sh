#!/bin/tcsh
setenv DATA /autofs/cluster/iaslab/FSMAP2/FSMAP_data
setenv SCRIPTPATH /autofs/cluster/iaslab/users/danlei/FSMAP/scripts
setenv IMAGE /autofs/cluster/iaslab/users/jtheriault/singularity_images/jtnipyutil/jtnipyutil-2019-01-03-4cecb89cb1d9.simg
setenv PROJNAME wm_nback_block
setenv SINGULARITY /usr/bin/singularity
setenv OUTPUT $DATA/BIDS_modeled/$PROJNAME

mkdir -p $OUTPUT/SC_warp
mkdir -p /scratch/$USER/$PROJNAME/wrkdir/
mkdir -p /scratch/$USER/$PROJNAME/output/
mkdir -p /scratch/$USER/$PROJNAME/data

rsync $SCRIPTPATH/model/res4d_SC_wrap/{wm_SC_mask_warp_startup.sh,wm_SC_mask_warp.py} /scratch/$USER/wm_nback_block/wrkdir/
chmod a+rwx /scratch/$USER/wm_nback_block/wrkdir/wm_SC_mask_warp_startup.sh
cd $home

$SINGULARITY exec  \
--bind "$DATA/BIDS_modeled/$PROJNAME/nosmooth:/scratch/data" \
--bind "$OUTPUT/SC_warp:/scratch/output" \
--bind "/scratch/$USER/$PROJNAME/wrkdir:/scratch/wrkdir" \
--bind "/autofs/cluster/iaslab/FSMAP2/FSMAP_data/BIDS_modeled/subject_SC_mask/wm/dartel_flow:/scratch/dartel_flow" \
$IMAGE \
/scratch/wrkdir/wm_SC_mask_warp_startup.sh

gzip $OUTPUT/SC_warp/*/*/*.nii

set ITERATIONS = (1 2 3 4 5 6 7 8 9 10 11 12)
foreach ITE ($ITERATIONS)
	echo ########################################################################
	echo ########################################################################
	echo $ITE
	echo ########################################################################
	echo ########################################################################
	##########################################
	$SINGULARITY exec  \
	--bind "$DATA/BIDS_modeled/$PROJNAME/nosmooth:/scratch/data" \
	--bind "$OUTPUT/SC_warp:/scratch/output" \
	--bind "/scratch/$USER/$PROJNAME/wrkdir:/scratch/wrkdir" \
	--bind "/autofs/cluster/iaslab/FSMAP2/FSMAP_data/BIDS_modeled/subject_SC_mask/wm/dartel_flow:/scratch/dartel_flow" \
	$IMAGE \
	/scratch/wrkdir/wm_SC_mask_warp_startup.sh
	gzip $OUTPUT/SC_warp/*/*/*.nii
	##########################################
	echo "start sleep"
	# sleep 18000 #18000=5hrs
	sleep 43200 
	#43200=12hours #18000=5hrs
	echo "end sleep"
end


# mv $OUTPUT/SC_warp/7_one_back_b1_cope/warped_files/w_sub-069__modelestimate0_cope4.nii.gz $OUTPUT/SC_warp/7_one_back_b1_cope/warped_files/w_sub-069__modelestimate0_cope7.nii.gz

# mkdir $OUTPUT
# rsync -r /scratch/$USER/$PROJNAME/output/ $OUTPUT

# cd /autofs/cluster/iaslab/users/danlei/FSMAP/scripts/model/
# chmod -R a+rwx *

rm -r /scratch/$USER/$PROJNAME/
exit


# scp -r * /autofs/cluster/iaslab/users/danlei/test/
# scp -r dz609@door.nmr.mgh.harvard.edu:/autofs/cluster/iaslab/users/danlei/test /Users/chendanlei/Desktop/
