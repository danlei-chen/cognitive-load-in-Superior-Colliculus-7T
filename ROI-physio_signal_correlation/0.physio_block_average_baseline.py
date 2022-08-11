import nipype.pipeline.engine as pe # pypeline engine
from nipype import IdentityInterface
import numpy as np
import pandas as pd
import os, sys
import nibabel as nib
import glob
import math

# subj_list = [i.split('/')[-1] for i in glob.glob('/Users/chendanlei/Desktop/U01/level1_files/wm_3vs1/SC_warp/*')]
# subj_list.sort()
subj_list = ['sub-012', 'sub-013', 'sub-015', 'sub-019', 'sub-021', 'sub-022', 'sub-023', 'sub-024', 'sub-025', 'sub-026', 'sub-030', 'sub-032', 'sub-034', 'sub-037', 'sub-039', 'sub-042b', 'sub-043', 'sub-044', 'sub-045', 'sub-047', 'sub-048', 'sub-049', 'sub-050', 'sub-052', 'sub-053', 'sub-054', 'sub-055', 'sub-056', 'sub-057', 'sub-058', 'sub-059', 'sub-060', 'sub-061', 'sub-062', 'sub-063', 'sub-064', 'sub-065', 'sub-066', 'sub-067', 'sub-069', 'sub-070', 'sub-071', 'sub-072', 'sub-073', 'sub-074', 'sub-078', 'sub-080', 'sub-081', 'sub-082', 'sub-083', 'sub-084', 'sub-085', 'sub-086', 'sub-087', 'sub-088', 'sub-090', 'sub-091', 'sub-092', 'sub-093', 'sub-094', 'sub-095', 'sub-098', 'sub-099', 'sub-100', 'sub-101', 'sub-102', 'sub-104', 'sub-105', 'sub-106', 'sub-110', 'sub-111', 'sub-112', 'sub-114', 'sub-117', 'sub-118', 'sub-119', 'sub-120b', 'sub-122b', 'sub-124', 'sub-127', 'sub-128', 'sub-131', 'sub-133', 'sub-134', 'sub-135', 'sub-137', 'sub-138', 'sub-139b']
event_files = glob.glob('/Volumes/GoogleDrive/My Drive/U01/working_memory/data/event_files/')
confound_files = glob.glob('/Volumes/GoogleDrive/My Drive/U01/working_memory/data/confound_files/')
save_dir = '/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/'

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

# physio_average_df = pd.DataFrame(columns=['block','block_type','block_type_number','block_total_number','block_TR_start','block_TR_end','subject'])
physio_average_df = pd.DataFrame()
for subj in subj_list:
    if subj[-1]=='b':
        subj = subj[0:7]
    print(subj)
    
    subj_average_df = pd.DataFrame(columns=['block','block_type','block_type_number','block_total_number','block_TR_start','block_TR_end','subject'])
    subj_df_event = pd.read_csv(glob.glob('/Volumes/GoogleDrive/My Drive/U01/working_memory/data/event_files/'+subj+'*.tsv')[0], sep='\t', header=0)
    subj_df_confound = pd.read_csv(glob.glob('/Volumes/GoogleDrive/My Drive/U01/working_memory/data/confound_files/'+subj+'*.tsv')[0], sep='\t', header=0)

    num_TR = subj_df_confound.shape[0]
    TR_series = [x*TR for x in range(0, num_TR)]  
    subj_df_confound['cumulative_TR'] = TR_series

    #block period
    subj_average_df['block_trial_start'] = [0,10,20,30,40,50,60,70,80,90,100,110]
    subj_average_df['block_trial_end'] = [i+9 for i in subj_average_df['block_trial_start']]
    
    block_TR_start_idx = [closest(TR_series, subj_df_event['onset'][i]) for i in subj_average_df['block_trial_start']]
    block_TR_end_idx = [closest(TR_series, subj_df_event['onset'][i]+subj_df_event['duration'][i]) for i in subj_average_df['block_trial_end']]
    subj_average_df['block_TR_start'] = [TR_series[i] for i in block_TR_start_idx]
    subj_average_df['block_TR_end'] = [TR_series[i] for i in block_TR_end_idx]
    
    #get baseline period
    baseline_TR_start_idx = [0]+[i+1 for i in block_TR_end_idx]
    baseline_TR_end_idx = [i-1 for i in block_TR_start_idx]+[num_TR]
    
    #log block type and other ingo    
    subj_average_df['block_type'] = [subj_df_event['trial_type'][i] for i in subj_average_df['block_trial_start']]
    subj_average_df.loc[subj_average_df['block_type']=='three_back','block_type_number'] = range(1,7)
    subj_average_df.loc[subj_average_df['block_type']=='one_back','block_type_number'] = range(1,7)
    subj_average_df['block'] = [i+'_b'+str(int(j)) for i, j in zip(subj_average_df['block_type'], subj_average_df['block_type_number'])]
    subj_average_df['block_total_number'] = range(1,len(subj_average_df['block_trial_start'])+1)
    subj_average_df['subject'] = subj
 
    try:
        EDA_DDA_nSCR_sum_baseline = np.nanmean([np.nansum(subj_df_confound['EDA_DDA_nSCR'][baseline_TR_start_idx[i] : baseline_TR_end_idx[i]]) for i in range(len(baseline_TR_start_idx))])
        EDA_DDA_tonic_mean_baseline = np.nanmean([np.nanmean(subj_df_confound['EDA_DDA_tonic'][baseline_TR_start_idx[i] : baseline_TR_end_idx[i]]) for i in range(len(baseline_TR_start_idx))])
        EDA_DDA_latency_mean_baseline = np.nanmean([np.nanmean(subj_df_confound['EDA_DDA_latency'][baseline_TR_start_idx[i] : baseline_TR_end_idx[i]]) for i in range(len(baseline_TR_start_idx))])
        EDA_DDA_ampsum_mean_baseline = np.nanmean([np.nanmean(subj_df_confound['EDA_DDA_ampsum'][baseline_TR_start_idx[i] : baseline_TR_end_idx[i]]) for i in range(len(baseline_TR_start_idx))])
        EDA_DDA_areasum_mean_baseline = np.nanmean([np.nanmean(subj_df_confound['EDA_DDA_areasum'][baseline_TR_start_idx[i] : baseline_TR_end_idx[i]]) for i in range(len(baseline_TR_start_idx))])
        EDA_global_mean_mean_baseline = np.nanmean([np.nanmean(subj_df_confound['EDA_global_mean'][baseline_TR_start_idx[i] : baseline_TR_end_idx[i]]) for i in range(len(baseline_TR_start_idx))])
        
        subj_average_df['EDA_DDA_nSCR_sum'] = [np.nansum(subj_df_confound['EDA_DDA_nSCR'][block_TR_start_idx[i] : block_TR_end_idx[i]]) for i in range(len(block_TR_start_idx))]
        subj_average_df['EDA_DDA_tonic_mean'] = [np.nanmean(subj_df_confound['EDA_DDA_tonic'][block_TR_start_idx[i] : block_TR_end_idx[i]]) for i in range(len(block_TR_start_idx))]
        subj_average_df['EDA_DDA_latency_mean'] = [np.nanmean(subj_df_confound['EDA_DDA_latency'][block_TR_start_idx[i] : block_TR_end_idx[i]]) for i in range(len(block_TR_start_idx))]
        subj_average_df['EDA_DDA_ampsum_mean'] = [np.nanmean(subj_df_confound['EDA_DDA_ampsum'][block_TR_start_idx[i] : block_TR_end_idx[i]]) for i in range(len(block_TR_start_idx))]
        subj_average_df['EDA_DDA_areasum_mean'] = [np.nanmean(subj_df_confound['EDA_DDA_areasum'][block_TR_start_idx[i] : block_TR_end_idx[i]]) for i in range(len(block_TR_start_idx))]
        subj_average_df['EDA_global_mean_mean'] = [np.nanmean(subj_df_confound['EDA_global_mean'][block_TR_start_idx[i] : block_TR_end_idx[i]]) for i in range(len(block_TR_start_idx))]

        subj_average_df['EDA_DDA_nSCR_sum'] =  subj_average_df['EDA_DDA_nSCR_sum'] - EDA_DDA_nSCR_sum_baseline
        subj_average_df['EDA_DDA_tonic_mean'] = subj_average_df['EDA_DDA_tonic_mean'] - EDA_DDA_tonic_mean_baseline
        subj_average_df['EDA_DDA_latency_mean'] = subj_average_df['EDA_DDA_latency_mean'] - EDA_DDA_latency_mean_baseline
        subj_average_df['EDA_DDA_ampsum_mean'] = subj_average_df['EDA_DDA_ampsum_mean'] - EDA_DDA_ampsum_mean_baseline
        subj_average_df['EDA_DDA_areasum_mean'] = subj_average_df['EDA_DDA_areasum_mean'] - EDA_DDA_areasum_mean_baseline
        subj_average_df['EDA_global_mean_mean'] = subj_average_df['EDA_global_mean_mean'] - EDA_global_mean_mean_baseline
        
    except:
        print('*********** no EDA ***********')
        
    try:
        resp_rate_smooth_mean_baseline = np.nanmean([np.nanmean(subj_df_confound['resp_rate_smooth'][baseline_TR_start_idx[i] : baseline_TR_end_idx[i]]) for i in range(len(baseline_TR_start_idx))])
        rvt_mean_baseline = np.nanmean([np.nanmean(subj_df_confound['rvt'][baseline_TR_start_idx[i] : baseline_TR_end_idx[i]]) for i in range(len(baseline_TR_start_idx))])

        subj_average_df['resp_rate_smooth_mean'] = [np.nanmean(subj_df_confound['resp_rate_smooth'][block_TR_start_idx[i] : block_TR_end_idx[i]]) for i in range(len(block_TR_start_idx))]
        subj_average_df['rvt_mean'] = [np.nanmean(subj_df_confound['rvt'][block_TR_start_idx[i] : block_TR_end_idx[i]]) for i in range(len(block_TR_start_idx))]
        
        subj_average_df['resp_rate_smooth_mean'] = subj_average_df['resp_rate_smooth_mean'] - resp_rate_smooth_mean_baseline
        subj_average_df['rvt_mean'] = subj_average_df['rvt_mean'] - rvt_mean_baseline
                
    except:
        print('*********** no respiration ***********')   

    try:
        ibi_mean_baseline = np.nanmean([np.nanmean(subj_df_confound['ibi'][baseline_TR_start_idx[i] : baseline_TR_end_idx[i]]) for i in range(len(baseline_TR_start_idx))])
        h_rate_mean_baseline = np.nanmean([np.nanmean(subj_df_confound['h_rate'][baseline_TR_start_idx[i] : baseline_TR_end_idx[i]]) for i in range(len(baseline_TR_start_idx))])

        subj_average_df['ibi_mean'] = [np.nanmean(subj_df_confound['ibi'][block_TR_start_idx[i] : block_TR_end_idx[i]]) for i in range(len(block_TR_start_idx))]
        subj_average_df['h_rate_mean'] = [np.nanmean(subj_df_confound['h_rate'][block_TR_start_idx[i] : block_TR_end_idx[i]]) for i in range(len(block_TR_start_idx))]

        subj_average_df['ibi_mean'] = subj_average_df['ibi_mean'] - ibi_mean_baseline
        subj_average_df['h_rate_mean'] = subj_average_df['h_rate_mean'] - h_rate_mean_baseline
        
    except:
        print('*********** no cardiac ***********')   
        
    physio_average_df = physio_average_df.append(subj_average_df,ignore_index=True)

physio_average_df.to_csv(save_dir+'physio_average_baseline.csv',index=True)
