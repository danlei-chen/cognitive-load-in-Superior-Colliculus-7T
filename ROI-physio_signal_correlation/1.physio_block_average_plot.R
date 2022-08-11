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
# df <- read_csv('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/physio_average.csv')
df <- read_csv('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/physio_average_baseline.csv')
# df <- read_csv('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/physio_average_baseline_percchange.csv')
fmri_subj <- unique(df$subject)
print('number of subjects: ')
print(length(fmri_subj))

#physio QA files
hr_qa <- read_csv('/Users/chendanlei/Dropbox (Partners HealthCare)/NCI_U01_Shared_Data/QA_documents/csv/FSMAP_valid_HR.csv')
hr_clean_subj <- hr_qa$subj[hr_qa$wm=='clean']
resp_qa <- read_csv('/Users/chendanlei/Dropbox (Partners HealthCare)/NCI_U01_Shared_Data/QA_documents/csv/FSMAP_valid_resp.csv')
resp_clean_subj <- resp_qa$subj[resp_qa$wm=='clean']
eda_qa <- read_csv('/Users/chendanlei/Dropbox (Partners HealthCare)/NCI_U01_Shared_Data/QA_documents/csv/FSMAP_valid_eda.csv')
eda_clean_subj <- eda_qa$subj[eda_qa$wm=='clean' | eda_qa$wm=='resp_noise']
eda_clean_subj <- eda_qa$subj[eda_qa$wm=='clean']
all_clean_subj <- intersect(intersect(hr_clean_subj, resp_clean_subj), eda_clean_subj)
  
######################################################################
#IBI, heart rate, RVT
######################################################################
ggplot(df[is.element(df$subject,hr_clean_subj),], aes(block_type_number, ibi_mean, colour = block_type)) + 
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
  scale_colour_manual(values = c('three_back' = "gray40", 'one_back' = "gray70")) + 
  theme(legend.position="none") +
  ggtitle(paste0('IBI, N=', length(unique(df[is.element(df$subject,hr_clean_subj),]$subject))))+ theme(text=element_text(family="Times New Roman"))
  # coord_cartesian(ylim = c(0, 6)) + 
anova(lmer(ibi_mean ~
                (block_type * block_type_number) +
                (1 | subject), data=df[is.element(df$subject,hr_clean_subj),]))
# summary(aov(ibi_mean~(block_type*block_type_number)+Error(factor(subject)), data = df))
t.test(df$ibi_mean[df$block_type=='three_back'], df$ibi_mean[df$block_type=='one_back'], paired = TRUE, alternative = "less")

ggplot(df[is.element(df$subject,hr_clean_subj),], aes(block_type_number, h_rate_mean, colour = block_type)) + 
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
  scale_colour_manual(values = c('three_back' = "gray40", 'one_back' = "gray70")) + 
  theme(legend.position="none") +
  ggtitle(paste0('heart rate, N=', length(unique(df[is.element(df$subject,hr_clean_subj),]$subject))))+ theme(text=element_text(family="Times New Roman"))
# coord_cartesian(ylim = c(0, 6)) + 
anova(lmer(h_rate_mean ~
             (block_type * block_type_number) +
             (1 | subject), data=df[is.element(df$subject,hr_clean_subj),]))
# summary(aov(h_rate_mean~(block_type*block_type_number)+Error(factor(subject)), data = df))
t.test(df$h_rate_mean[df$block_type=='three_back'], df$h_rate_mean[df$block_type=='one_back'], paired = TRUE, alternative = "less")

######################################################################
#respiration
######################################################################
ggplot(df[is.element(df$subject,resp_clean_subj),], aes(block_type_number, resp_rate_smooth_mean, colour = block_type)) + 
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
  scale_colour_manual(values = c('three_back' = "gray40", 'one_back' = "gray70")) + 
  theme(legend.position="none") +
  ggtitle(paste0('resp rate, N=', length(unique(df[is.element(df$subject,resp_clean_subj),]$subject))))+ theme(text=element_text(family="Times New Roman"))
# coord_cartesian(ylim = c(0, 6)) + 
anova(lmer(resp_rate_smooth_mean ~
                (block_type * block_type_number) +
                (1 | subject), data=df[is.element(df$subject,hr_clean_subj),]))
# summary(aov(rvt_mean~(block_type*block_type_number)+Error(factor(subject)), data = df))
t.test(df$resp_rate_smooth_mean[df$block_type=='three_back'], df$resp_rate_smooth_mean[df$block_type=='one_back'], paired = TRUE, alternative = "greater")

ggplot(df[is.element(df$subject,resp_clean_subj),], aes(block_type_number, rvt_mean, colour = block_type)) + 
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
  scale_colour_manual(values = c('three_back' = "gray40", 'one_back' = "gray70")) + 
  theme(legend.position="none") +
  ggtitle(paste0('block mean RVT, N=', length(unique(df[is.element(df$subject,resp_clean_subj),]$subject))))+ theme(text=element_text(family="Times New Roman"))
# coord_cartesian(ylim = c(0, 6)) + 
anova(lmer(rvt_mean ~
             (block_type + block_type_number) +
             (1 | subject), data=df[is.element(df$subject,resp_clean_subj),]))
# summary(aov(rvt_mean~(block_type*block_type_number)+Error(factor(subject)), data = df))
t.test(df$rvt_mean[df$block_type=='three_back'], df$rvt_mean[df$block_type=='one_back'], paired = TRUE, alternative = "greater")

######################################################################
#EDA
######################################################################
########
# event
ggplot(df[is.element(df$subject,eda_clean_subj),], aes(block_type_number, EDA_DDA_nSCR_sum, colour = block_type)) + 
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
  scale_colour_manual(values = c('three_back' = "gray40", 'one_back' = "gray70")) + 
  theme(legend.position="none") +
  ggtitle(paste0('SCR events, N=', length(unique(df[is.element(df$subject,eda_clean_subj),]$subject))))+ theme(text=element_text(family="Times New Roman"))
ggplot(df[is.element(df$subject,eda_clean_subj),], aes(block_type, EDA_DDA_nSCR_sum, fill = block_type)) + 
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
  scale_colour_manual(values = c('three_back' = "gray40", 'one_back' = "gray70")) + 
  ggtitle('block mean SCR events')+ theme(text=element_text(family="Times New Roman"))
# coord_cartesian(ylim = c(0, 6)) + 
# theme(legend.position="none")
anova(lmer(EDA_DDA_nSCR_sum ~
                (block_type + block_type_number) +
                (1 | subject), data=df[is.element(df$subject,eda_clean_subj),]))  
# summary(aov(rvt_mean~(block_type*block_type_number)+Error(factor(subject)), data = df))
t.test(df$EDA_DDA_nSCR_sum[df$block_type=='three_back'], df$EDA_DDA_nSCR_sum[df$block_type=='one_back'], paired = TRUE, alternative = "greater")

ggplot(df[is.element(df$subject,eda_clean_subj),], aes(block_type_number, EDA_DDA_tonic_mean, colour = block_type)) + 
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
  scale_colour_manual(values = c('three_back' = "gray40", 'one_back' = "gray70")) + 
  theme(legend.position="none") +
  ggtitle(paste0('tonic EDA mean, N=', length(unique(df[is.element(df$subject,eda_clean_subj),]$subject))))+ theme(text=element_text(family="Times New Roman"))
# coord_cartesian(ylim = c(0, 6)) + 
anova(lmer(EDA_DDA_tonic_mean ~
                (block_type * block_type_number) +
                (1 | subject), data=df[is.element(df$subject,eda_clean_subj),]))  
# summary(aov(rvt_mean~(block_type*block_type_number)+Error(factor(subject)), data = df))
t.test(df$EDA_DDA_tonic_mean[df$block_type=='three_back'], df$EDA_DDA_tonic_mean[df$block_type=='one_back'], paired = TRUE, alternative = "greater")

ggplot(df[is.element(df$subject,eda_clean_subj),], aes(block_type_number, EDA_DDA_ampsum_mean, colour = block_type)) + 
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
  scale_colour_manual(values = c('three_back' = "gray40", 'one_back' = "gray70")) + 
  theme(legend.position="none") +
  ggtitle(paste0('SCR amplitude sum, N=', length(unique(df[is.element(df$subject,eda_clean_subj),]$subject))))+ theme(text=element_text(family="Times New Roman"))
# coord_cartesian(ylim = c(0, 6)) + 
# theme(legend.position="none")
anova(lmer(EDA_DDA_ampsum_mean ~
             (block_type + block_type_number) +
             (1 | subject), data=df[is.element(df$subject,eda_clean_subj),]))  
# summary(aov(rvt_mean~(block_type*block_type_number)+Error(factor(subject)), data = df))
t.test(df$EDA_DDA_ampsum_mean[df$block_type=='three_back'], df$EDA_DDA_ampsum_mean[df$block_type=='one_back'], paired = TRUE, alternative = "greater")

ggplot(df[is.element(df$subject,eda_clean_subj),], aes(block_type_number, EDA_DDA_areasum_mean, colour = block_type)) + 
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
  scale_colour_manual(values = c('three_back' = "gray40", 'one_back' = "gray70")) + 
  theme(legend.position="none") +
  ggtitle(paste0('SCR area sum, N=', length(unique(df[is.element(df$subject,eda_clean_subj),]$subject))))+ theme(text=element_text(family="Times New Roman"))
# coord_cartesian(ylim = c(0, 6)) + 
# theme(legend.position="none")
anova(lmer(EDA_DDA_areasum_mean ~
             (block_type + block_type_number) +
             (1 | subject), data=df[is.element(df$subject,eda_clean_subj),]))  
# summary(aov(rvt_mean~(block_type*block_type_number)+Error(factor(subject)), data = df))
t.test(df$EDA_DDA_areasum_mean[df$block_type=='three_back'], df$EDA_DDA_areasum_mean[df$block_type=='one_back'], paired = TRUE, alternative = "greater")

ggplot(df[is.element(df$subject,eda_clean_subj),], aes(block_type_number, EDA_global_mean_mean, colour = block_type)) + 
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
  scale_colour_manual(values = c('three_back' = "gray40", 'one_back' = "gray70")) + 
  theme(legend.position="none") +
  ggtitle(paste0('EDA global mean, N=', length(unique(df[is.element(df$subject,eda_clean_subj),]$subject))))+ theme(text=element_text(family="Times New Roman"))
# coord_cartesian(ylim = c(0, 6)) + 
# theme(legend.position="none")
anova(lmer(EDA_global_mean_mean ~
             (block_type + block_type_number) +
             (1 | subject), data=df[is.element(df$subject,eda_clean_subj),]))  
# summary(aov(rvt_mean~(block_type*block_type_number)+Error(factor(subject)), data = df))
t.test(df$EDA_global_mean_mean[df$block_type=='three_back'], df$EDA_global_mean_mean[df$block_type=='one_back'], paired = TRUE, alternative = "greater")

ggplot(df[is.element(df$subject,eda_clean_subj),], aes(block_type_number, EDA_DDA_latency_mean, colour = block_type)) + 
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
  scale_colour_manual(values = c('three_back' = "gray40", 'one_back' = "gray70")) + 
  theme(legend.position="none") +
  ggtitle('block mean DDA latency')+ theme(text=element_text(family="Times New Roman"))
# coord_cartesian(ylim = c(0, 6)) + 
# theme(legend.position="none")
anova(lmer(EDA_DDA_latency_mean ~
             (block_type * block_type_number) +
             (block_type * block_type_number | subject), data=df[is.element(df$subject,eda_clean_subj),]))  
# summary(aov(rvt_mean~(block_type*block_type_number)+Error(factor(subject)), data = df))
t.test(df$EDA_DDA_latency_mean[df$block_type=='three_back'], df$EDA_DDA_latency_mean[df$block_type=='one_back'], paired = TRUE, alternative = "greater")

