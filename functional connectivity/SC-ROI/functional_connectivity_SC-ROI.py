#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 13:10:46 2021

@author: chendanlei
"""

#python3 /Volumes/GoogleDrive/My\ Drive/U01/working_memory/analysis/connectivity/scripts/functional_connectivity_SC-ROI.py

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
smooth_mm_list = [3]
roi_dir_base = '/scratch/data/ROI_mask/'
rois = [os.environ['ROI']]
subj_sc_dir = '/scratch/data/SC_mask/'
hemispheres = ['all', 'rh', 'lh']
template = '/scratch/data/res4d/w_sub-012__modelestimate0_res4d.nii.gz'

output_dir_base = '/scratch/output/'
# os.makedirs(output_dir_base, exist_ok = True)
input_dir = '/scratch/data/res4d/'
subj_list = [i.split('w_')[1].split('__')[0] for i in glob.glob(input_dir+'*')]
# subj_list=['sub-042','sub-012']
subj_list.sort()
print(subj_list)
print(len(subj_list))

confound_file_dir = '/scratch/data/confound_files/'
event_file_dir = '/scratch/data/confound_files/'
#find the index of the closest number to K in a list
def closest(lst, K): 
    val = lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))] 
    idx = lst.index(val)
    return idx
#get consecutive intervals of 
def get_consecutive_intervals(nums):
    nums = sorted(set(nums))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s+1 < e]
    edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
    return list(zip(edges, edges))
TR=2.34

roi_df = pd.DataFrame(columns=['subject','seed','target','smooth','hemisphere','block_type','corr_coef'])

for roi in rois:
    print('SC to '+roi)
    
    for hemi in hemispheres:
        print('hemispheres: '+hemi)
        
        #get and resample roi mask
        roi_img_orig = nib.load(glob.glob(roi_dir_base+'*'+roi+'_'+hemi+'.nii.gz')[0])
        roi_img = resample_img(roi_img_orig,
                                target_affine=nib.load(template).affine,
                                target_shape=nib.load(template).shape[0:3],
                                interpolation='nearest')
        
        for smooth in smooth_mm_list:
            print('smoothing : '+str(smooth)+' mm')
            
            for subj in subj_list:
                print(subj)
                
                #get, resample, and threshold subject sc mask
                subj_sc_img_orig = nib.load(glob.glob(subj_sc_dir+'*'+subj+'*')[0])
                subj_sc_data = resample_img(subj_sc_img_orig,
                                        target_affine=nib.load(template).affine,
                                        target_shape=nib.load(template).shape[0:3],
                                        interpolation='nearest').get_fdata()
                subj_sc_data[subj_sc_data<.25] = 0
                subj_sc_data[subj_sc_data!=0] = 1
                if hemi=='lh':
                    subj_sc_data[range(int(subj_sc_data.shape[0]/2), int(subj_sc_data.shape[0])),:,:]=0
                elif hemi=='rh':
                    subj_sc_data[range(0, int(subj_sc_data.shape[0]/2)),:,:]=0
                subj_sc_img = nib.Nifti1Image(subj_sc_data, nib.load(template).affine, nib.load(template).header)

                #get res4d file 
                res4d_img = nib.load(glob.glob(input_dir+'/w_'+subj+'*res4d.nii.gz')[0])
                
                #smooth res4d
                res4d_img_smoothed = smooth_img(res4d_img,fwhm=smooth)
            
                #mask res 4d sc
                sc_masker = NiftiMasker(mask_img=subj_sc_img, standardize=False)
                sc_masker.fit(res4d_img_smoothed) 
                res4d_img_smoothed_sc_masked_mean = np.mean(sc_masker.transform(res4d_img_smoothed),axis=1)
        
                #mask res 4d roi
                roi_masker = NiftiMasker(mask_img=roi_img, standardize=False)
                roi_masker.fit(res4d_img_smoothed) 
                res4d_img_smoothed_roi_masked_mean = np.mean(roi_masker.transform(res4d_img_smoothed),axis=1)
                
                ############ correlation of all blocks ############
                #find correlation 
                corr_coef_all = np.corrcoef(res4d_img_smoothed_sc_masked_mean, res4d_img_smoothed_roi_masked_mean)[0,1]
                roi_df = roi_df.append({'subject':subj,'seed':'SC','target':roi,'smooth':smooth,'hemisphere':hemi,'block_type':'all','corr_coef':corr_coef_all}, ignore_index=True)

                ############ correlation of 3-/1-back blocks ############
                #get block condition 
                subj_df_event = pd.read_csv(glob.glob(event_file_dir+subj+'*/func/'+subj+'*_task-wm_events.tsv')[0], sep='\t', header=0)
                subj_df_confound = pd.read_csv(glob.glob(confound_file_dir+subj+'*/func/'+subj+'*_task-wm_bold_confounds.tsv')[0], sep='\t', header=0)
            
                num_TR = subj_df_confound.shape[0]
                TR_series = [x*TR for x in range(0, num_TR)]  
            
                #block period
                block_trial_start = [0,10,20,30,40,50,60,70,80,90,100,110]
                block_trial_end = [i+9 for i in block_trial_start]
                # block_trial_start_end = [(0,9),(10,19),(20,29),(30,39),(40,49),(50,59),(60,69),(70,79),(80,89),(90,99),(100,109),(110,119)]

                block_TR_start_idx = [closest(TR_series, subj_df_event['onset'][i]) for i in block_trial_start]
                block_TR_end_idx = [closest(TR_series, subj_df_event['onset'][i]+subj_df_event['duration'][i]) for i in block_trial_end]
                block_TR_idx = [(block_TR_start_idx[i], block_TR_end_idx[i]) for i in range(0, len(block_TR_start_idx))]

                block_condition = [subj_df_event['trial_type'][i] for i in block_trial_start]
                
                res4d_img_smoothed_sc_masked_mean_3b = [list(res4d_img_smoothed_sc_masked_mean[block_TR_idx[i][0] : block_TR_idx[i][1]]) for i in [i for i, x in enumerate(block_condition) if x == 'three_back']]
                res4d_img_smoothed_sc_masked_mean_3b = [item for sublist in res4d_img_smoothed_sc_masked_mean_3b for item in sublist]
                res4d_img_smoothed_roi_masked_mean_3b = [list(res4d_img_smoothed_roi_masked_mean[block_TR_idx[i][0] : block_TR_idx[i][1]]) for i in [i for i, x in enumerate(block_condition) if x == 'three_back']]
                res4d_img_smoothed_roi_masked_mean_3b = [item for sublist in res4d_img_smoothed_roi_masked_mean_3b for item in sublist]
                #find correlation 
                corr_coef_3b = np.corrcoef(res4d_img_smoothed_sc_masked_mean_3b, res4d_img_smoothed_roi_masked_mean_3b)[0,1]
                roi_df = roi_df.append({'subject':subj,'seed':'SC','target':roi,'smooth':smooth,'hemisphere':hemi,'block_type':'3-back','corr_coef':corr_coef_3b}, ignore_index=True)

                res4d_img_smoothed_sc_masked_mean_1b = [list(res4d_img_smoothed_sc_masked_mean[block_TR_idx[i][0] : block_TR_idx[i][1]]) for i in [i for i, x in enumerate(block_condition) if x == 'one_back']]
                res4d_img_smoothed_sc_masked_mean_1b = [item for sublist in res4d_img_smoothed_sc_masked_mean_1b for item in sublist]
                res4d_img_smoothed_roi_masked_mean_1b = [list(res4d_img_smoothed_roi_masked_mean[block_TR_idx[i][0] : block_TR_idx[i][1]]) for i in [i for i, x in enumerate(block_condition) if x == 'one_back']]
                res4d_img_smoothed_roi_masked_mean_1b = [item for sublist in res4d_img_smoothed_roi_masked_mean_1b for item in sublist]
                #find correlation 
                corr_coef_1b = np.corrcoef(res4d_img_smoothed_sc_masked_mean_1b, res4d_img_smoothed_roi_masked_mean_1b)[0,1]
                roi_df = roi_df.append({'subject':subj,'seed':'SC','target':roi,'smooth':smooth,'hemisphere':hemi,'block_type':'1-back','corr_coef':corr_coef_1b}, ignore_index=True)
                
                roi_df.to_csv(output_dir_base+'/SC-'+roi+'_subj.csv')







                # corr_map_img = nib.Nifti1Image(corr_map, nib.load(con).affine, nib.load(con).header)
                # corr_map_img_name = conn_output_dir+cope_name+'_'+str(smooth_mm)+'mm_'+conn_map_name_suffix+'_corrcoef'+seed_suffix+'.nii.gz'
                # nib.save(corr_map_img, corr_map_img_name)
    
                # corr_map_fisher_transform = np.arctanh(corr_map)
                # corr_map_fisher_transform_img = nib.Nifti1Image(corr_map_fisher_transform, nib.load(con).affine, nib.load(con).header)
                # corr_map_fisher_transform_img_name = conn_output_dir+cope_name+'_'+str(smooth_mm)+'mm_'+conn_map_name_suffix+'_corrcoef_FisherZ'+seed_suffix+'.nii.gz'
                # nib.save(corr_map_fisher_transform_img, corr_map_fisher_transform_img_name)

                
    
