#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 13:10:46 2021

@author: chendanlei
"""

#python3 /Volumes/GoogleDrive/My\ Drive/U01/working_memory/analysis/connectivity/scripts/SCvoxel-ROI/SCvoxel_average.py

import os
import nibabel as nib
import numpy as np
import pandas as pd
import glob
# from nilearn.glm import threshold_stats_img
from nilearn.image import resample_img
from nilearn.image import smooth_img
from nilearn.input_data import NiftiMasker

smooth_mm_list = [None, 1.5, 3]
rois = ['FEF', 'LIP', 'LGN', 'V1', 'MGN', 'A1', 'LIPd', 'LIPv']
sc_dir = '/Volumes/GoogleDrive/My Drive/U01/working_memory/roi/subject_SC_mask/avg_template/Template_6_thresh0.25_resammpled_bin.nii.gz'
input_dir_base = '/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/connectivity/results/SCvoxel-ROI/'
output_dir_base = '/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/connectivity/results/SCvoxel-ROI/'
hemispheres = ['all', 'rh', 'lh']
template = '/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/connectivity/data/res4d_SC_warp/wm_nback_block/sub-012/w_sub-012__modelestimate0_res4d.nii.gz'
distance_info = pd.read_csv('/Volumes/GoogleDrive/My Drive/U01/working_memory/roi/SC_layer_by_distance/sc_segment_info_cope_space.csv', index_col=0)

template_data = nib.load(template).dataobj[:,:,:,0]
#empty csv
roi_mean_df = pd.DataFrame(columns=['seed','seed_xyz','seed_x','seed_y','seed_z','target','smooth','hemisphere','block_type','corr_coef_Z_mean', 'distance', 'side', 'quartile'])
                
for roi in rois:
    print('SC voxel to '+roi)
    
    roi_df = pd.read_csv(input_dir_base+'SCvoxel-'+roi+'_subj.csv', index_col=0)
    roi_df['corr_coef_Z'] = np.arctanh(np.array(roi_df['corr_coef']))
    
    for hemi in hemispheres:
        print('hemispheres: '+hemi)

        for smooth in smooth_mm_list:
            print('smoothing : '+str(smooth)+' mm')
            
            for condition in ['3-back','1-back','all']:
            
                #empty template
                average_data = template_data.copy()
                average_data[:]=0

                if smooth==None:
                    roi_df_sub = roi_df[(roi_df['hemisphere']==hemi) & (np.isnan(roi_df['smooth'])) & (roi_df['block_type']==condition)]
                else:
                    roi_df_sub = roi_df[(roi_df['hemisphere']==hemi) & (roi_df['smooth']==smooth) & (roi_df['block_type']==condition)]
                    
                for seed_xyz in np.unique(roi_df_sub['seed_xyz']):
                    x = int(seed_xyz[1:3])
                    y = int(seed_xyz[5:7])
                    z = int(seed_xyz[9:11])
                    
                    roi_df_sub_vox = roi_df_sub[(roi_df_sub['seed_x']==x) & (roi_df_sub['seed_y']==y) & (roi_df_sub['seed_z']==z)]
            
                    average_data[x, y, z] = np.mean(roi_df_sub_vox['corr_coef_Z'])
                    roi_mean_df = roi_mean_df.append({'seed':'SC','seed_xyz':(x,y,z),'seed_x':x,'seed_y':y,'seed_z':z,'target':roi,'smooth':smooth,'hemisphere':hemi,'block_type':condition,'corr_coef_Z_mean':average_data[x, y, z],'distance':float(distance_info['distance'][(distance_info['x']==x) & (distance_info['y']==y) & (distance_info['z']==z)]), 'side':list(distance_info['side'][(distance_info['x']==x) & (distance_info['y']==y) & (distance_info['z']==z)])[0], 'quartile':int(distance_info['quartile'][(distance_info['x']==x) & (distance_info['y']==y) & (distance_info['z']==z)])}, ignore_index=True)

                nib.save(nib.Nifti1Image(average_data, nib.load(template).affine, nib.load(template).header), os.path.join(output_dir_base, 'average', 'SCvoxel-'+roi+'_'+hemi+'_'+str(smooth)+'mm_mean.nii.gz'))
                roi_mean_df.to_csv(output_dir_base+'/SCvoxel-ROI_mean.csv')
