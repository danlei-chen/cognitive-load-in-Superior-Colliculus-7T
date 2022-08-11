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

df_cor <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(df_cor) <- c("subject","subregion", "type",  "physio", "roi", "physio_signal", "roi_signal")

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
    print(physio_variable)
    for (subregion in c('all')){
      for (subj in subject_intersect){
        # print(subj)
        for (type in c("3-back","1-back")){
          roi_sig <- as.double(df_roi$mean_signal[df_roi$subject==subj & df_roi$type==type & df_roi$roi_dim1==subregion])
          physio_sig <- as.double(as.matrix(df_physio[physio_variable])[df_physio$subject==subj & df_physio$block_type==type])
          if (subj %in% c("sub-069", "sub-080")){
            physio_sig <- physio_sig[1:3]
          }
          for (n in c(1:length(physio_sig))){
            df_cor[nrow(df_cor)+1,] <- c(as.character(subj), as.character(subregion), as.character(type), as.character(physio_variable), as.character(roi), as.double(physio_sig[n]), as.double(roi_sig[n]))
          }
        }
      }
    }
  }
}
write.csv(df_cor, '/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/physio-roi_average_correlation_baseline2.csv')

df_cor <- read_csv('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/physio-roi_average_correlation_baseline2.csv')
# df_cor <- read_csv('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/physio-roi_average_correlation_baseline_percchange2.csv')
df_cor$roi_signal <- as.double(df_cor$roi_signal)
df_cor$physio_signal <- as.double(df_cor$physio_signal)
# df_cor$roi_signal_zscore <- (df_cor$roi_signal-mean(df_cor$roi_signal))/sd(df_cor$roi_signal)
# df_cor$physio_signal_zscore <- (df_cor$physio_signal-mean(df_cor$physio_signal, na.rm = TRUE))/sd(df_cor$physio_signal, na.rm = TRUE)

df_cor$roi <- factor(df_cor$roi, levels=c('SC','SN','LGN','Caudate','LPBN','MPBN','VSM'))
df_cor$physio = factor(df_cor$physio, levels=c("IBI","heart_rate","resp_rate","RVT","nSCR","tonic_EDA","EDA_latency","EDA_ampsum","EDA_areasum","EDA_global_mean"))

#mean subject correlation
# roi_list = c('SC', 'SN', 'Caudate', 'MPBN')
# roi_list = c('SC', 'SN', 'Caudate', 'LPBN', 'MPBN', 'VSM', 'LGN')
roi_list = c('SC')
# df_cor <- df_cor[df_cor$physio == 'IBI' | df_cor$physio == 'resp_rate' | df_cor$physio == 'EDA_ampsum',]

#physio QA files
hr_qa <- read_csv('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/physio_QA_csv/FSMAP_valid_HR.csv')
hr_clean_subj <- hr_qa$subj[hr_qa$wm=='clean']
resp_qa <- read_csv('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/physio_QA_csv/FSMAP_valid_resp.csv')
resp_clean_subj <- resp_qa$subj[resp_qa$wm=='clean']
eda_qa <- read_csv('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/physio_QA_csv/FSMAP_valid_eda.csv')
eda_clean_subj <- eda_qa$subj[eda_qa$wm=='clean' | eda_qa$wm=='resp_noise']
# eda_clean_subj <- eda_qa$subj[eda_qa$wm=='clean']
all_clean_subj <- intersect(intersect(hr_clean_subj, resp_clean_subj), eda_clean_subj)

df_cor_clean <- rbind(
  df_cor[df_cor$physio=='IBI',][df_cor$subject[df_cor$physio=='IBI'] %in% hr_clean_subj,], 
  df_cor[df_cor$physio=='resp_rate',][df_cor$subject[df_cor$physio=='resp_rate'] %in% resp_clean_subj,], 
  df_cor[df_cor$physio=='EDA_ampsum',][df_cor$subject[df_cor$physio=='EDA_ampsum'] %in% eda_clean_subj,],
  df_cor[df_cor$physio=='tonic_EDA',][df_cor$subject[df_cor$physio=='tonic_EDA'] %in% eda_clean_subj,])

df_sub <- df_cor_clean[df_cor_clean$physio=='IBI' & df_cor_clean$roi=='SC' & df_cor_clean$subregion=='all',]
# df_sub <- df_cor_clean[df_cor_clean$type=='3-back' & df_cor_clean$physio=='IBI' & df_cor_clean$roi=='SC' & df_cor_clean$subregion=='all',]
ggplot(df_sub, aes(roi_signal, physio_signal, colour = type)) + 
  stat_summary(fun = mean, geom = "point", size=1, alpha=0.5) +
  # stat_summary(fun = mean, geom = "line", size=.25, alpha=0.5) +
  geom_smooth(method=lm)
anova(lmer(physio_signal ~ (roi_signal * type) + (roi_signal * type | subject), data=df_sub))
anova(lmer(roi_signal ~ (physio_signal * type) + (physio_signal * type | subject), data=df_sub))

df_sub <- df_cor_clean[df_cor_clean$physio=='resp_rate' & df_cor_clean$roi=='SC' & df_cor_clean$subregion=='all',]
ggplot(df_sub, aes(roi_signal, physio_signal, colour = type)) + 
  stat_summary(fun = mean, geom = "point", size=1, alpha=0.5) +
  # stat_summary(fun = mean, geom = "line", size=.25, alpha=0.5) +
  geom_smooth(method=lm)
anova(lmer(physio_signal ~ (roi_signal * type) + (roi_signal * type | subject), data=df_sub))
anova(lmer(roi_signal ~ (physio_signal * type) + (physio_signal * type | subject), data=df_sub))

df_sub <- df_cor_clean[df_cor_clean$physio=='EDA_ampsum' & df_cor_clean$roi=='SC' & df_cor_clean$subregion=='all',]
tonic_EDA_sig <- df_cor_clean$physio_signal[df_cor_clean$physio=='tonic_EDA' & df_cor_clean$roi=='SC' & df_cor_clean$subregion=='all']
df_sub <- cbind(df_sub, tonic_EDA_sig)
ggplot(df_sub, aes(roi_signal, physio_signal, colour = type)) + 
  stat_summary(fun = mean, geom = "point", size=1, alpha=0.5) +
  # stat_summary(fun = mean, geom = "line", size=.25, alpha=0.5) +
  geom_smooth(method=lm)
anova(lmer(roi_signal ~
             ((physio_signal*tonic_EDA_sig) * type | subject), data=df_sub))
anova(lmer(physio_signal ~
             ((roi_signal) * type) +
             ((roi_signal) * type | subject), data=df_sub))


#mediation analysis
# https://towardsdatascience.com/doing-and-reporting-your-first-mediation-analysis-in-r-2fe423b92171
df_sub <- df_cor_clean[df_cor_clean$physio=='EDA_ampsum' & df_cor_clean$roi=='SC' & df_cor_clean$subregion=='all',]
fit.totaleffect=lm(roi_signal~type,df_sub)
summary(fit.totaleffect)
fit.mediator=lm(physio_signal~type,df_sub)
summary(fit.mediator)
fit.dv=lm(roi_signal~physio_signal+type,df_sub)
summary(fit.dv)

library(mediation)
results = mediate(fit.mediator, fit.dv, treat='type', mediator='physio_signal', boot=T)
summary(results)








for (roi in roi_list){
  print(roi)
  for (physio in unique(df_cor_clean$physio)){
    print(physio)
    
    plot <- ggplot(df_cor_clean[df_cor_clean$roi==roi & df_cor_clean$physio==physio,], aes(x = roi_signal, y = physio_signal, group=subject, colour = subject)) +
    geom_point(alpha=0.3)+
    geom_smooth(method = "lm", se = FALSE, alpha=0.6, size=.5)+
    facet_grid(~type)+
    ggtitle(paste0(roi,' ',physio))
    # theme(legend.position="none")
    print(plot)
    
    # plot <- ggplot(df_cor[df_cor$roi==roi & df_cor$physio==physio,], aes(x = roi_signal_zscore, y = physio_signal_zscore, group=subject, colour = subject)) +
    #   geom_point(alpha=0.3)+
    #   geom_smooth(method = "lm", se = FALSE, alpha=0.6, size=.5)+
    #   facet_grid(~type)+
    #   ggtitle(paste0(roi,' ',physio))
    # # theme(legend.position="none")
    # print(plot)
    
  }
}

