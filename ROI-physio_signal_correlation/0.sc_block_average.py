#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 21:25:00 2021

@author: chendanlei
"""
#chmod a+rwx /Volumes/GoogleDrive/My\ Drive/U01/working_memory/analysis/physio/0.sc_block_average.py 
#python3 /Volumes/GoogleDrive/My\ Drive/U01/working_memory/analysis/physio/0.sc_block_average.py 

import nibabel as nib
import numpy as np
import glob
from nilearn.glm import threshold_stats_img
from nilearn.image import resample_img
import pandas as pd 
import os
import statistics

#subject level roi directory
roi='SC'
data_dir = '/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/univariate/results/level1/wm_nback_block/SC_warp/'
roi_dir = '/Volumes/GoogleDrive/My Drive/U01/working_memory/roi/subject_SC_mask/warped_SC/'
roi_thresh = 0.25
output = '/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/'+roi+'_signal_wm_'+str(roi_thresh)+'_qua.csv'
task='wm'
# file_format1 = 'cope2_smoothed3mm_masked.nii.gz'
# file_format2 = 'cope3_smoothed3mm_masked.nii.gz'
file_format_list = ['cope'+str(i)+'.nii.gz' for i in range(1,13)]
file_format_list_name = ['3-back_b'+str(i) for i in range(1,7)]+['1-back_b'+str(i) for i in range(1,7)]
template_file = '/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/univariate/results/level1/wm_nback_block/SC_warp/sub-012/w_sub-012__modelestimate0_cope1.nii.gz'

subj_list = [i.split('/')[-1] for i in glob.glob(data_dir+'*')]
subj_list = [i[0:7] for i in subj_list if 'b' in i]
subj_list.sort()

df = pd.DataFrame(columns = ['subjID', 'subject','task','type','block','file_name','roi_dim1','roi_dim3','roi_size','mean_signal','peak_signal'])
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
        
            #load roi
            # roi_file = glob.glob(roi_dir+'/*SC*'+s+'*')[0]
            roi_file = glob.glob(roi_dir+'*'+s+'*')[0]
            roi_data = nib.load(roi_file).get_fdata()    
            # roi_data = resample_img(nib.load(roi_file),
            #     target_affine=file_img.affine,
            #     target_shape=file_img.shape[0:3],
            #     interpolation='nearest').get_fdata()    
            roi_data[roi_data<roi_thresh] = np.nan
            roi_data[pd.notnull(roi_data)] = 1
            roi_size = len(roi_data[pd.notnull(roi_data)])
            
            #mask contrast file
            file_data[np.isnan(roi_data)] = np.nan
            #get mean
            mean_signal = np.nanmean(file_data)
            #get peak
            max_signal = np.nanmax(file_data)
            #write to df
            df_row = pd.Series([s, subj, task, file_format_list_name[file_format_num].split('_b')[0], file_format_list_name[file_format_num].split('_b')[-1], file_format_list_name[file_format_num], 'all', 'all', roi_size,mean_signal, max_signal], index = df.columns)
            df = df.append(df_row, ignore_index=True)
            
            #get upper left roi
            roi_data = nib.load(roi_file).get_fdata()    
            # roi_data = resample_img(nib.load(roi_file),
            #     target_affine=file_img.affine,
            #     target_shape=file_img.shape[0:3],
            #     interpolation='nearest').get_fdata()    
            roi_data[roi_data<roi_thresh] = np.nan
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
            df_row = pd.Series([s, subj, task, file_format_list_name[file_format_num].split('_b')[0], file_format_list_name[file_format_num].split('_b')[-1], file_format_list_name[file_format_num], 'left', 'upper', roi_size, mean_signal, max_signal], index = df.columns)
            df = df.append(df_row, ignore_index=True)
        
            #get lower left roi
            roi_data = nib.load(roi_file).get_fdata()    
            # roi_data = resample_img(nib.load(roi_file),
            #     target_affine=file_img.affine,
            #     target_shape=file_img.shape[0:3],
            #     interpolation='nearest').get_fdata()    
            roi_data[roi_data<roi_thresh] = np.nan
            roi_data[pd.notnull(roi_data)] = 1
            roi_data[int(roi_data.shape[0]/2):,:,:] = np.nan #first dim second half - left,
            height_middle_split = int(statistics.median(np.where(~np.isnan(roi_data))[2]))
            roi_data[:,:,height_middle_split:] = np.nan #third dim second half - lower,
            roi_size = len(roi_data[pd.notnull(roi_data)])
            #mask contrast file
            file_data = nib.load(file).get_fdata()
            file_data[np.isnan(roi_data)] = np.nan
            #get mean
            mean_signal = np.nanmean(file_data)
            #get peak
            max_signal = np.nanmax(file_data)
            #write to df
            df_row = pd.Series([s, subj, task, file_format_list_name[file_format_num].split('_b')[0], file_format_list_name[file_format_num].split('_b')[-1], file_format_list_name[file_format_num], 'left', 'lower', roi_size, mean_signal, max_signal], index = df.columns)
            df = df.append(df_row, ignore_index=True)
        
            #get upper right roi
            roi_data = nib.load(roi_file).get_fdata()    
            # roi_data = resample_img(nib.load(roi_file),
            #     target_affine=file_img.affine,
            #     target_shape=file_img.shape[0:3],
            #     interpolation='nearest').get_fdata()    
            roi_data[roi_data<roi_thresh] = np.nan
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
            df_row = pd.Series([s, subj, task, file_format_list_name[file_format_num].split('_b')[0], file_format_list_name[file_format_num].split('_b')[-1], file_format_list_name[file_format_num], 'right', 'upper', roi_size, mean_signal, max_signal], index = df.columns)
            df = df.append(df_row, ignore_index=True)
            
            #get lower right roi
            roi_data = nib.load(roi_file).get_fdata()    
            # roi_data = resample_img(nib.load(roi_file),
            #     target_affine=file_img.affine,
            #     target_shape=file_img.shape[0:3],
            #     interpolation='nearest').get_fdata()    
            roi_data[roi_data<roi_thresh] = np.nan
            roi_data[pd.notnull(roi_data)] = 1
            roi_data[0:int(roi_data.shape[0]/2),:,:] = np.nan #first dim second half - right
            height_middle_split = int(statistics.median(np.where(~np.isnan(roi_data))[2]))
            roi_data[:,:,height_middle_split:] = np.nan #third dim second half - lower,
            roi_size = len(roi_data[pd.notnull(roi_data)])
            #mask contrast file
            file_data = nib.load(file).get_fdata()
            file_data[np.isnan(roi_data)] = np.nan
            #get mean
            mean_signal = np.nanmean(file_data)
            #get peak
            max_signal = np.nanmax(file_data)
            #write to df
            df_row = pd.Series([s, subj, task, file_format_list_name[file_format_num].split('_b')[0], file_format_list_name[file_format_num].split('_b')[-1], file_format_list_name[file_format_num], 'right', 'lower', roi_size, mean_signal, max_signal], index = df.columns)
            df = df.append(df_row, ignore_index=True)

        except:
            print(subj+' has a problem in '+file_format)
    
    df.to_csv(output)
    
    
    