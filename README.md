# cognitive-load-in-Superior-Colliculus-7T

The project focuses on modeling superior colliculus responses during a working memory task with 2 cognitive loads (3-back and 1-back). Data was collected in high-resolution 7 Tesla fMRI (N=88).

The scripts include subject-level extration of superior colliculus masks;, univariate analysis of task modeling and warping of first level result on group-level superior colliculus mask, MVPA analysis using SVM to classify cognitive loads, and roi and physiological signal extraction and correlation analysis.

Main analyses were done in Python (Sklearn, Nipype, Nilearn), mixed effects statistical model was done in R (lmer), visualization was mostly done in R (ggplot).
