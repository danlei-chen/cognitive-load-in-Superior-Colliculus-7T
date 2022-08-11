#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 21:25:00 2021

@author: chendanlei
"""
#chmod a+rwx /Volumes/GoogleDrive/My\ Drive/U01/working_memory/analysis/physio/0.roi_block_average.py 
#python3 /Volumes/GoogleDrive/My\ Drive/U01/working_memory/analysis/physio/0.roi_block_average.py 

import nibabel as nib
import numpy as np
import glob
from nilearn.glm import threshold_stats_img
from nilearn.image import resample_img
import pandas as pd 
import os
import statistics

#subject level roi directory
roi_list=['SN', 'Caudate', 'LPBN', 'MPBN', 'VSM', 'LGN']
data_dir = '/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/univariate/results/level1/wm_nback_block/SC_warp/'
task='wm'
file_format_list = ['cope'+str(i)+'.nii.gz' for i in range(1,13)]
file_format_list_name = ['3-back_b'+str(i) for i in range(1,7)]+['1-back_b'+str(i) for i in range(1,7)]
template_file = '/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/univariate/results/level1/wm_nback_block/SC_warp/sub-012/w_sub-012__modelestimate0_cope1.nii.gz'

subj_list = [i.split('/')[-1] for i in glob.glob(data_dir+'*')]
subj_list = [i[0:7] for i in subj_list if 'b' in i]
subj_list.sort()

for roi in roi_list:
    print(roi)
    output = '/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/'+roi+'_signal_wm.csv'
    
    if roi=='LGN':
        roi_dir = '/Volumes/GoogleDrive/My Drive/fMRI/roi/7Trep_Marta/LG_l_r_resammpled.nii.gz'
    elif roi=='SN':
        roi_dir = '/Volumes/GoogleDrive/My Drive/fMRI/roi/7Trep_Marta/SN_l_r_resammpled.nii.gz'
    elif roi=='LPBN':
        roi_dir = '/Volumes/GoogleDrive/My Drive/fMRI/roi/7Trep_Marta/LPB_l_r_resammpled.nii.gz'
    elif roi=='MPBN':
        roi_dir = '/Volumes/GoogleDrive/My Drive/fMRI/roi/7Trep_Marta/MPB_l_r_resammpled.nii.gz'
    elif roi=='VSM':
        roi_dir = '/Volumes/GoogleDrive/My Drive/fMRI/roi/7Trep_Marta/VSM_l_r_resammpled.nii.gz'
    elif roi=='Caudate':
        roi_dir = '/Volumes/GoogleDrive/My Drive/fMRI/roi/Pauli_thresh05_combined/Caudate.nii.gz'
        
    #load roi
    roi_data_orig = resample_img(nib.load(roi_dir),
        target_affine=nib.load(template_file).affine,
        target_shape=nib.load(template_file).shape[0:3],
        interpolation='nearest').get_fdata()    
        
    df = pd.DataFrame(columns = ['subjID', 'subject','task','type','block','file_name','roi_dim1','roi_size','mean_signal','peak_signal'])
    for s in subj_list:
        print(s)
        
        subj = s
        
        for file_format_num,file_format in enumerate(file_format_list):

            try:

                ######## NEGATIVE COPE #######s
                #load contrast
                file = glob.glob(data_dir+s+'*/*'+file_format)[0]
                file_img = nib.load(file)
                file_data = file_img.get_fdata()
            
                # #load roi
                # roi_data = resample_img(nib.load(roi_dir),
                #     target_affine=file_img.affine,
                #     target_shape=file_img.shape[0:3],
                #     interpolation='nearest').get_fdata() 
                roi_data = roi_data_orig.copy()
                roi_data[roi_data!=1] = np.nan
                roi_data[pd.notnull(roi_data)] = 1
                roi_size = len(roi_data[pd.notnull(roi_data)])
                
                #mask contrast file
                file_data[np.isnan(roi_data)] = np.nan
                #get mean
                mean_signal = np.nanmean(file_data)
                #get peak
                max_signal = np.nanmax(file_data)
                #write to df
                df_row = pd.Series([s, subj, task, file_format_list_name[file_format_num].split('_b')[0], file_format_list_name[file_format_num].split('_b')[-1], file_format_list_name[file_format_num], 'all', roi_size,mean_signal, max_signal], index = df.columns)
                df = df.append(df_row, ignore_index=True)
                
                # #get left roi
                # roi_data = resample_img(nib.load(roi_dir),
                #     target_affine=file_img.affine,
                #     target_shape=file_img.shape[0:3],
                #     interpolation='nearest').get_fdata()  
                roi_data = roi_data_orig.copy()
                roi_data[roi_data!=1] = np.nan
                roi_data[pd.notnull(roi_data)] = 1
                roi_data[int(roi_data.shape[0]/2):,:,:] = np.nan #first dim second half - left,
                height_middle_split = int(statistics.median(np.where(~np.isnan(roi_data))[2]))
                roi_data[:,:,0:height_middle_split] = np.nan #third dim first half - upper,
                roi_size = len(roi_data[pd.notnull(roi_data)])
                #mask contrast file
                file_data = nib.load(file).get_fdata()
                file_data[np.isnan(roi_data)] = np.nan
                # test_img = nib.Nifti1Image(file_data, nib.load(file).affine, nib.load(file).header)
                # nib.save(test_img, '/Users/chendanlei/Desktop/x.nii.gz')
                #get mean
                mean_signal = np.nanmean(file_data)
                #get peak
                max_signal = np.nanmax(file_data)
                #write to df
                df_row = pd.Series([s, subj, task, file_format_list_name[file_format_num].split('_b')[0], file_format_list_name[file_format_num].split('_b')[-1], file_format_list_name[file_format_num], 'left', roi_size, mean_signal, max_signal], index = df.columns)
                df = df.append(df_row, ignore_index=True)
            
                # #get right roi
                # roi_data = resample_img(nib.load(roi_dir),
                #     target_affine=file_img.affine,
                #     target_shape=file_img.shape[0:3],
                #     interpolation='nearest').get_fdata()  
                roi_data = roi_data_orig.copy()
                roi_data[roi_data!=1] = np.nan
                roi_data[pd.notnull(roi_data)] = 1
                roi_data[0:int(roi_data.shape[0]/2),:,:] = np.nan #first dim second half - right
                height_middle_split = int(statistics.median(np.where(~np.isnan(roi_data))[2]))
                roi_data[:,:,0:height_middle_split] = np.nan #third dim first half - upper,
                roi_size = len(roi_data[pd.notnull(roi_data)])
                #mask contrast file
                file_data = nib.load(file).get_fdata()
                file_data[np.isnan(roi_data)] = np.nan
                #get mean
                mean_signal = np.nanmean(file_data)
                #get peak
                max_signal = np.nanmax(file_data)
                #write to df
                df_row = pd.Series([s, subj, task, file_format_list_name[file_format_num].split('_b')[0], file_format_list_name[file_format_num].split('_b')[-1], file_format_list_name[file_format_num], 'right', roi_size, mean_signal, max_signal], index = df.columns)
                df = df.append(df_row, ignore_index=True)

                df.to_csv(output)
            
            except:
                print(subj+' has a problem in '+file_format)


        
    
        