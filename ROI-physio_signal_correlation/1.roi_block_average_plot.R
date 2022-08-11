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
# roi_list = c('SC', 'SN', 'Caudate', 'LPBN', 'MPBN', 'VSM', 'LGN')
roi_list = c('SN', 'Caudate', 'LPBN', 'MPBN', 'VSM', 'LGN')
for (roi in roi_list){
  print(roi)
  tryCatch(
    {
      df <- read_csv(paste0('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/',roi,'_signal_wm.csv'))
    },error=function(e){
      df <- read_csv(paste0('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/',roi,'_signal_wm_0.25_qua.csv'))
    }
  )
  print('number of subjects: ')
  print(length(unique(df$subject)))

  plot <- ggplot(df[df$roi_dim1=='all',], aes(block, mean_signal, colour = type)) + 
    # geom_point(alpha=0.2)+
    stat_summary(fun = mean, geom = "point") +
    stat_summary(fun = mean, geom = "line") +
    stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
    scale_colour_manual(values = c('3-back' = "gray40", '1-back' = "gray70")) + 
    ggtitle(paste0('all', ' ', roi, ' block mean signal'))+ theme(text=element_text(family="Times New Roman"))
  print(plot)
  # coord_cartesian(ylim = c(0, 6)) + 
  # theme(legend.position="none")
  print(anova(lmer(mean_signal ~ (type + block) + (type + block | subject), data=df[df$roi_dim1=='all',])))
  # summary(aov(ibi_mean~(block_type*block_type_number)+Error(factor(subject)), data = df))
  # print(t.test(df$mean_signal[df$type=='3-back' & df$roi_dim1=='all'], df$mean_signal[df$type=='1-back' & df$roi_dim1=='all'], paired = TRUE, alternative = "less"))
  
  plot <- ggplot(df[df$roi_dim1=='all',], aes(type, mean_signal, fill = type)) + 
    # geom_point(alpha=0.2)+
    stat_summary(fun = mean, geom = "bar") +
    stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
    scale_fill_manual(values = c('3-back' = "gray40", '1-back' = "gray70")) + 
    ggtitle(paste0(roi, ' block mean signal'))+ theme(text=element_text(family="Times New Roman"))
    # facet_wrap(~roi_dim1)
  print(plot)  
  print(t.test(df$mean_signal[df$type=='3-back' & df$roi_dim1=='all'], df$mean_signal[df$type=='1-back' & df$roi_dim1=='all'], paired = TRUE, alternative = "greater"))
  
  plot <- ggplot(df[df$roi_dim1!='all',], aes(block, mean_signal, colour = type)) + 
    # geom_point(alpha=0.2)+
    stat_summary(fun = mean, geom = "point") +
    stat_summary(fun = mean, geom = "line") +
    stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
    scale_colour_manual(values = c('3-back' = "gray40", '1-back' = "gray70")) + 
    ggtitle(paste0(roi, ' block mean signal'))+ theme(text=element_text(family="Times New Roman")) +
    facet_wrap(~roi_dim1)
  print(plot)
  # coord_cartesian(ylim = c(0, 6)) + 
  # theme(legend.position="none")
  print(anova(lmer(mean_signal ~ (type + block + roi_dim1) + (type + block + roi_dim1| subject), data=df[df$roi_dim1!='all',])))
  # summary(aov(ibi_mean~(block_type*block_type_number)+Error(factor(subject)), data = df))
  
  plot <- ggplot(df[df$roi_dim1!='all',], aes(type, mean_signal, fill = type)) + 
    # geom_point(alpha=0.2)+
    stat_summary(fun = mean, geom = "bar") +
    stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
    scale_fill_manual(values = c('3-back' = "gray40", '1-back' = "gray70")) + 
    ggtitle(paste0(roi, ' block mean signal'))+ theme(text=element_text(family="Times New Roman")) +
    facet_wrap(~roi_dim1)
  print(plot)  
  print(t.test(df$mean_signal[df$type=='3-back' & df$roi_dim1!='all'], df$mean_signal[df$type=='1-back' & df$roi_dim1!='all'], paired = TRUE, alternative = "greater"))
  
  rm(df)
  
}

######################################################################
######################################################################
roi = 'SC'
print(roi)
df <- read_csv(paste0('/Volumes/GoogleDrive/My Drive/U01/working_memory/analysis/physio/block_average_results/',roi,'_signal_wm_0.25_qua.csv'))
print('number of subjects: ')
print(length(unique(df$subject)))
df$roi_dim3 <- factor(df$roi_dim3, levels = c('upper', 'lower'))

plot <- ggplot(df[df$roi_dim1=='all',], aes(block, mean_signal, colour = type)) + 
  # geom_point(alpha=0.2)+
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
  scale_colour_manual(values = c('3-back' = "gray40", '1-back' = "gray70")) + 
  ggtitle(paste0('all', ' ', roi, ' block mean signal'))+ theme(text=element_text(family="Times New Roman"))
print(plot)
# coord_cartesian(ylim = c(0, 6)) + 
# theme(legend.position="none")
print(anova(lmer(mean_signal ~ (type + block) + (type + block | subject), data=df[df$roi_dim1=='all',])))
# summary(aov(ibi_mean~(block_type*block_type_number)+Error(factor(subject)), data = df))
# print(t.test(df$mean_signal[df$type=='3-back' & df$roi_dim1=='all'], df$mean_signal[df$type=='1-back' & df$roi_dim1=='all'], paired = TRUE, alternative = "less"))

plot <- ggplot(df[df$roi_dim1=='all',], aes(type, mean_signal, fill = type)) + 
  # geom_point(alpha=0.2)+
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
  scale_fill_manual(values = c('3-back' = "gray40", '1-back' = "gray70")) + 
  ggtitle(paste0(roi, ' block mean signal'))+ theme(text=element_text(family="Times New Roman"))
# facet_wrap(~roi_dim1)
print(plot)  
print(t.test(df$mean_signal[df$type=='3-back' & df$roi_dim1=='all'], df$mean_signal[df$type=='1-back' & df$roi_dim1=='all'], paired = TRUE, alternative = "greater"))

plot <- ggplot(df[df$roi_dim1!='all',], aes(block, mean_signal, colour = type)) + 
  # geom_point(alpha=0.2)+
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
  scale_colour_manual(values = c('3-back' = "gray40", '1-back' = "gray70")) + 
  ggtitle(paste0(roi, ' block mean signal'))+ theme(text=element_text(family="Times New Roman")) +
  facet_wrap(roi_dim3~roi_dim1)
print(plot)
# coord_cartesian(ylim = c(0, 6)) + 
# theme(legend.position="none")
print(anova(lmer(mean_signal ~ (type + block + roi_dim1 + roi_dim3) + (type + block + roi_dim1 + roi_dim3| subject), data=df[df$roi_dim1!='all',])))
# summary(aov(ibi_mean~(block_type*block_type_number)+Error(factor(subject)), data = df))

plot <- ggplot(df[df$roi_dim1!='all',], aes(type, mean_signal, fill = type)) + 
  # geom_point(alpha=0.2)+
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom="errorbar",width=0.2, size = .5) +
  scale_fill_manual(values = c('3-back' = "gray40", '1-back' = "gray70")) + 
  ggtitle(paste0(roi, ' block mean signal'))+ theme(text=element_text(family="Times New Roman")) +
  facet_wrap(roi_dim3~roi_dim1)
print(plot)  
print(t.test(df$mean_signal[df$type=='3-back' & df$roi_dim1!='all'], df$mean_signal[df$type=='1-back' & df$roi_dim1!='all'], paired = TRUE, alternative = "greater"))

