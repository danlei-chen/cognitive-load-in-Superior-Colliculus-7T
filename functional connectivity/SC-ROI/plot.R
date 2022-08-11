library(ggplot2)
library(readr)
library(tidyr)
library(lme4)
library(nlme)
library(DescTools)

rois <- c('FEF', 'LIP', 'LGN', 'V1', 'MGN', 'A1', 'LIPd', 'LIPv')
rois <- c('FEF', 'LIP', 'LGN', 'V1', 'MGN', 'A1')
df_all <- data.frame()
for (roi in rois){
  df <- read_csv(paste0('/Volumes/GoogleDrive/My\ Drive/U01/working_memory/analysis/connectivity/results/SC-ROI/SC-',roi,'_subj.csv'))
  df_all <- rbind(df_all, df)
}
df_all$target <- factor(df_all$target, levels=c('FEF', 'LIP', 'LGN', 'MGN', 'V1', 'A1'))
df_all$corr_coef_Z <- FisherZ(df_all$corr_coef)

df_all <- df_all[is.na(df_all$smooth),]
df_all <- df_all[df_all$hemisphere=='all',]

ggplot(df_all[df_all$block_type=='all' & df_all$hemisphere=='all',], aes(x=target,y=corr_coef_Z)) + 
  stat_summary(fun=mean,geom="bar", size = 0.5, alpha=0.9)+
  stat_summary(fun=mean, fun.data=mean_se, geom="errorbar",width=0.2,size = 0.5)+
  # geom_jitter(alpha = 0.8, size = 1.5, width = 0.05, height = 0.05)+
  geom_point(alpha = 0.1, size = 1.5, width = 0, height = 0)+
  geom_line(alpha = 0.1, aes(group = subject))+
  # geom_hline(aes(yintercept=0.5), linetype='dotted') + 
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=13),
        axis.title=element_text(size=15),
        plot.title = element_text(size=16,face="bold"))

ggplot(df_all[df_all$block_type=='3-back'|df_all$block_type=='1-back',], aes(target, corr_coef_Z, group=block_type, fill=block_type)) +
  # geom_point(alpha=0.2, aes(colour=subject))+
  stat_summary(fun = mean, geom = "col", position = position_dodge()) +
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.2, size = .5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c('3-back' = "gray40", '1-back' = "gray70")) +
  # geom_hline(yintercept=0, linetype='dashed', alpha=0.5)+
  # ggtitle(paste0(roi)) +
  # facet_grid(block_type~hemisphere~smooth)
  # facet_grid(block_type~smooth)
  facet_grid(block_type~hemisphere)

ggplot(df_all[df_all$block_type=='3-back'|df_all$block_type=='1-back',], aes(target, corr_coef_Z, group=block_type, fill=block_type)) +
  # geom_point(alpha=0.2, aes(colour=subject))+
  stat_summary(fun = mean, geom = "col", position = position_dodge()) +
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.2, size = .5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c('3-back' = "gray40", '1-back' = "gray70")) +
  # geom_hline(yintercept=0, linetype='dashed', alpha=0.5)+
  # ggtitle(paste0(roi)) +
  facet_grid(hemisphere~smooth)

t.test(df_all$accuracy[df_all$k_folds=='88 folds'], mu = 0.5, alternative = "greater")
t.test(df_all$accuracy[df_all$k_folds=='12 folds'], mu = 0.5, alternative = "greater")
t.test(df_all$accuracy[df_all$k_folds=='6 folds'], mu = 0.5, alternative = "greater")
mean(df_all$accuracy[df_all$k_folds=='88 folds'])
sd(df_all$accuracy[df_all$k_folds=='88 folds'])/sqrt(88)

