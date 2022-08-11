library(ggplot2)
library(readr)
library(tidyr)
library(lme4)
library(nlme)

######### task classification########
# df_4folds <- read_tsv('/Volumes/GoogleDrive/My\ Drive/U01/working_memory/analysis/MVPA/1vs3back/1vs3back_MVPA_classification_6folds.tsv')
# df_10folds <- read_tsv('/Volumes/GoogleDrive/My\ Drive/U01/working_memory/analysis/MVPA/1vs3back/1vs3back_MVPA_classification_12folds.tsv')
# df_80folds <- read_tsv('/Volumes/GoogleDrive/My\ Drive/U01/working_memory/analysis/MVPA/1vs3back/1vs3back_MVPA_classification_88folds.tsv')
# df_all <- data.frame(accuracy = c(df_4folds$`6 folds`, df_10folds$`12 folds`, df_80folds$`88 folds`), k_folds = c(rep('6 folds',nrow(df_4folds)),rep('12 folds',nrow(df_10folds)),rep('88 folds',nrow(df_80folds))))
# df_all$k_folds <- factor(df_all$k_folds, levels = c('6 folds', '12 folds', '88 folds'))
df_80folds <- read_tsv('/Volumes/GoogleDrive/My\ Drive/U01/working_memory/analysis/MVPA/1vs3back/1vs3back_MVPA_classification_88folds.tsv')
df_all <- data.frame(accuracy = c(df_80folds$`88 folds`), k_folds = c(rep('88 folds',nrow(df_80folds))))
df_all$k_folds <- factor(df_all$k_folds, levels = c('88 folds'))

ggplot(df_all, aes(x=k_folds,y=accuracy)) + 
  stat_summary(fun=mean, alpha = 0.6, geom="bar", size = 0.5)+
  stat_summary(fun=mean, fun.data=mean_se, geom="errorbar",width=0.2,size = 0.5)+
  # geom_jitter(alpha = 0.8, size = 1.5, width = 0.05, height = 0.05)+
  geom_jitter(alpha = 0.2, size = 1.5, width = 0.2, height = 0)+
  geom_hline(aes(yintercept=0.5), linetype='dotted') + 
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=13),
        axis.title=element_text(size=15),
        plot.title = element_text(size=16,face="bold"))

t.test(df_all$accuracy[df_all$k_folds=='88 folds'], mu = 0.5, alternative = "greater")
t.test(df_all$accuracy[df_all$k_folds=='12 folds'], mu = 0.5, alternative = "greater")
t.test(df_all$accuracy[df_all$k_folds=='6 folds'], mu = 0.5, alternative = "greater")
mean(df_all$accuracy[df_all$k_folds=='88 folds'])
sd(df_all$accuracy[df_all$k_folds=='88 folds'])/sqrt(88)




