library(ggplot2)
library(readr)
library(tidyr)
library(lme4)
library(nlme)
library(lmerTest)
library(pbkrtest)
library(emmeans)

######################################################################
######################################################################
# df_physio <- read_csv('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/physio_average.csv')
df_physio <- read_csv('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/physio_average_baseline.csv')
# df_physio <- read_csv('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/physio_average_baseline_percchange.csv')
print('number of subjects in physio: ')
print(length(unique(df_physio$subject)))
df_physio$block_type[df_physio$block_type=='three_back']='3-back'
df_physio$block_type[df_physio$block_type=='one_back']='1-back'
colnames(df_physio) <- c("X1","block","block_type","block_type_number","block_total_number","block_TR_start","block_TR_end","subject","block_trial_start","block_trial_end","nSCR","tonic_EDA","EDA_latency","EDA_ampsum","EDA_areasum","EDA_global_mean","resp_rate","IBI","RVT","heart_rate")

roi_list = c('SC', 'SN', 'Caudate', 'LPBN', 'MPBN', 'VSM', 'LGN')
physio_variable_list = c("IBI","heart_rate","resp_rate","RVT","nSCR","tonic_EDA","EDA_latency","EDA_ampsum","EDA_areasum","EDA_global_mean")

df_mean <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(df_mean) <- c("subject", "physio", "roi",  "subregion", "type", "ROI_mean", "physio_mean")

for (roi in roi_list){
  print(roi)
  if (file.exists(paste0('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/',roi,'_signal_wm.csv'))){
    df_roi <- read_csv(paste0('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/',roi,'_signal_wm.csv'))
  }else{
    df_roi <- read_csv(paste0('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/',roi,'_signal_wm_0.25_qua.csv'))
  }

  subject_intersect <- intersect(unique(df_roi$subject),unique(df_physio$subject))
  print('number of subjects in both: ')
  print(length(subject_intersect))
  
  df_physio <- df_physio[df_physio$subject %in% subject_intersect,]
  df_roi <- df_roi[df_roi$subject %in% subject_intersect,]

  for (physio_variable in physio_variable_list){
    for (subregion in c('all', 'left', 'right')){
      for (subj in subject_intersect){
        # print(subj)
        for (type in c("3-back","1-back")){
          roi_sig <- as.double(df_roi$mean_signal[df_roi$subject==subj & df_roi$type==type & df_roi$roi_dim1==subregion])
          physio_sig <- as.double(as.matrix(df_physio[physio_variable])[df_physio$subject==subj & df_physio$block_type==type])
          if (subj %in% c("sub-069", "sub-080")){
            physio_sig <- physio_sig[1:3]
          }
          # if (length(roi_sig)> length(physio_sig)){
          #   physio_sig <-rep(physio_sig, each=2)
          # }
          df_mean[nrow(df_mean)+1,] <- c(as.character(subj), as.character(physio_variable), as.character(roi), as.character(subregion), as.character(type), as.double(mean(roi_sig)), as.double(mean(physio_sig)))
        }
      }
    }
  }
}
write.csv(df_mean, '/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/physio-roi_average_baseline_mean.csv')

##############################################################################################################
#read in csv
##############################################################################################################
df_mean <- read_csv('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/physio-roi_average_baseline_mean.csv')

#physio QA files
hr_qa <- read_csv('/Users/chendanlei/Dropbox (Partners HealthCare)/NCI_U01_Shared_Data/QA_documents/csv/FSMAP_valid_HR.csv')
hr_clean_subj <- hr_qa$subj[hr_qa$wm=='clean']
resp_qa <- read_csv('/Users/chendanlei/Dropbox (Partners HealthCare)/NCI_U01_Shared_Data/QA_documents/csv/FSMAP_valid_resp.csv')
resp_clean_subj <- resp_qa$subj[resp_qa$wm=='clean']
eda_qa <- read_csv('/Users/chendanlei/Dropbox (Partners HealthCare)/NCI_U01_Shared_Data/QA_documents/csv/FSMAP_valid_eda.csv')
eda_clean_subj <- eda_qa$subj[eda_qa$wm=='clean' | eda_qa$wm=='resp_noise']
# eda_clean_subj <- eda_qa$subj[eda_qa$wm=='clean']
all_clean_subj <- intersect(intersect(hr_clean_subj, resp_clean_subj), eda_clean_subj)

df_mean_clean <- rbind(df_mean[df_mean$physio=='IBI',][df_mean$subject[df_mean$physio=='IBI'] %in% hr_clean_subj,], df_mean[df_mean$physio=='resp_rate',][df_mean$subject[df_mean$physio=='resp_rate'] %in% resp_clean_subj,], df_mean[df_mean$physio=='EDA_ampsum',][df_mean$subject[df_mean$physio=='EDA_ampsum'] %in% eda_clean_subj,])

ggplot(df_mean_clean[df_mean_clean$roi=="SC",], aes(ROI_mean, physio_mean, color=type)) +
  geom_point(alpha=0.2)+
  geom_smooth(method=lm)+
  facet_grid(physio~subregion,scales="free")

##############################################################################################################
#cross subject correlation
##############################################################################################################
df_cor <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df_cor) <- c("physio", "roi",  "subregion", "type", "correlation")
for (roi in unique(df_mean_clean$roi)){
  for (physio_variable in unique(df_mean_clean$physio)){
    for (subregion in unique(df_mean_clean$subregion)){
      for (type in unique(df_mean_clean$type)){
        cor_value <- cor(df_mean_clean$physio_mean[df_mean_clean$roi==roi & df_mean_clean$physio==physio_variable & df_mean_clean$type==type & df_mean_clean$subregion==subregion],
                         df_mean_clean$ROI_mean[df_mean_clean$roi==roi & df_mean_clean$physio==physio_variable & df_mean_clean$type==type & df_mean_clean$subregion==subregion], 
            use = "pairwise.complete.obs")
        df_cor[nrow(df_cor)+1,] <- c(as.character(physio_variable), as.character(roi), as.character(subregion), as.character(type), as.double(cor_value))
}}}}
df_cor$correlation <- as.numeric(df_cor$correlation)
# # roi_list = c('SC', 'SN', 'Caudate', 'MPBN')
roi_list = c('SN', 'Caudate', 'LPBN', 'MPBN', 'VSM', 'LGN', 'SC')
df_cor <- df_cor[df_cor$physio == 'IBI' | df_cor$physio == 'resp_rate' | df_cor$physio == 'EDA_ampsum',]

ggplot(df_cor[df_cor$roi=='SC',], aes(type, correlation, fill=type)) +
  # geom_point(alpha=0.2, aes(colour=subject))+
  stat_summary(fun.y = mean, geom = "bar") +
  # stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.1, size = .5) +
  scale_fill_manual(values = c('3-back' = "gray40", '1-back' = "gray70")) + 
  geom_hline(yintercept=0, linetype='dashed', alpha=0.5)+
  ggtitle(paste0('cross-subject correlation between physiological signals and SC')) +
  facet_grid(physio~subregion)

ggplot(df_cor, aes(type, correlation, fill=type)) +
    # geom_point(alpha=0.2, aes(colour=subject))+
    stat_summary(fun.y = mean, geom = "bar") +
    # stat_summary(fun.y = mean, geom = "line") +
    stat_summary(fun.data = mean_se, geom="errorbar", width=0.1, size = .5) +
    scale_fill_manual(values = c('3-back' = "gray40", '1-back' = "gray70")) + 
    geom_hline(yintercept=0, linetype='dashed', alpha=0.5)+
    ggtitle(paste0('cross-subject correlation between physiological signals and ', roi)) +
    facet_grid(physio~subregion~roi)

ggplot(df_cor, aes(roi, correlation, group=1, fill=roi)) +
  # geom_point(alpha=0.2, aes(colour=subject))+
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.1, size = .5) +
  geom_hline(yintercept=0, linetype='dashed', alpha=0.5)+
  # ggtitle(paste0(roi)) +
  facet_grid(physio~subregion)

##############################################################################################################
#find 3-back and 1-back difference
##############################################################################################################
df_mean_clean_diff <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df_mean_clean_diff) <- c("subject", "physio", "roi",  "subregion", "ROI_mean_difference", "physio_mean_difference")
for(subj in unique(df_mean_clean$subject)){
  for (roi in unique(df_mean_clean$roi)){
    for (physio_variable in unique(df_mean_clean$physio)){
      for (subregion in unique(df_mean_clean$subregion)){
        physio_diff <- df_mean_clean$physio_mean[df_mean_clean$type=="3-back" & df_mean_clean$subject==subj & df_mean_clean$roi==roi & df_mean_clean$physio==physio_variable & df_mean_clean$subregion==subregion] -
          df_mean_clean$physio_mean[df_mean_clean$type=="1-back" & df_mean_clean$subject==subj & df_mean_clean$roi==roi & df_mean_clean$physio==physio_variable & df_mean_clean$subregion==subregion]
        roi_diff <- df_mean_clean$ROI_mean[df_mean_clean$type=="3-back" & df_mean_clean$subject==subj & df_mean_clean$roi==roi & df_mean_clean$physio==physio_variable & df_mean_clean$subregion==subregion] -
          df_mean_clean$ROI_mean[df_mean_clean$type=="1-back" & df_mean_clean$subject==subj & df_mean_clean$roi==roi & df_mean_clean$physio==physio_variable & df_mean_clean$subregion==subregion]
        if (length(physio_diff)==1){
          df_mean_clean_diff[nrow(df_mean_clean_diff)+1,] <- c(as.character(subj), as.character(physio_variable), as.character(roi), as.character(subregion), as.double(roi_diff), as.double(physio_diff))
        }
      }}}}
df_mean_clean_diff$ROI_mean_difference <- as.numeric(df_mean_clean_diff$ROI_mean_difference)
df_mean_clean_diff$physio_mean_difference <- as.numeric(df_mean_clean_diff$physio_mean_difference)
ggplot(df_mean_clean_diff[df_mean_clean_diff$roi=="SC",], aes(ROI_mean_difference, physio_mean_difference)) +
  geom_point(alpha=0.2)+
  geom_smooth(method=lm)+
  facet_grid(physio~subregion,scales="free")

df_cor <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(df_cor) <- c("physio", "roi",  "subregion", "correlation")
for (roi in unique(df_mean_clean_diff$roi)){
  for (physio_variable in unique(df_mean_clean_diff$physio)){
    for (subregion in unique(df_mean_clean_diff$subregion)){
      cor_value <- cor(df_mean_clean_diff$physio_mean[df_mean_clean_diff$roi==roi & df_mean_clean_diff$physio==physio_variable & df_mean_clean_diff$subregion==subregion],
                       df_mean_clean_diff$ROI_mean[df_mean_clean_diff$roi==roi & df_mean_clean_diff$physio==physio_variable & df_mean_clean_diff$subregion==subregion], 
                       use = "pairwise.complete.obs")
      df_cor[nrow(df_cor)+1,] <- c(as.character(physio_variable), as.character(roi), as.character(subregion), as.double(cor_value))
      }}}
df_cor$correlation <- as.numeric(df_cor$correlation)

ggplot(df_cor[df_cor$roi=="SC",], aes(subregion, correlation)) +
  # geom_point(alpha=0.2, aes(colour=subject))+
  stat_summary(fun.y = mean, geom = "bar") +
  # stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.1, size = .5) +
  scale_fill_manual(values = c('3-back' = "gray40", '1-back' = "gray70")) + 
  geom_hline(yintercept=0, linetype='dashed', alpha=0.5)+
  ggtitle(paste0('cross-subject correlation between physiological signals and ', roi)) +
  facet_grid(physio~1)

ggplot(df_cor, aes(subregion, correlation)) +
  # geom_point(alpha=0.2, aes(colour=subject))+
  stat_summary(fun.y = mean, geom = "bar") +
  # stat_summary(fun.y = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.1, size = .5) +
  scale_fill_manual(values = c('3-back' = "gray40", '1-back' = "gray70")) + 
  geom_hline(yintercept=0, linetype='dashed', alpha=0.5)+
  ggtitle(paste0('cross-subject correlation between physiological signals and ', roi)) +
  facet_grid(physio~roi)
