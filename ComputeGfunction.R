###This code estimates the G-functions for a specified COCO CDM and WDM sample
###The estimated G-functions used in our paper can be loaded directly from 
####./data/coco_mw_functional_summaries/
###Permutation p-values can be computed using ComputePermutationPvalues.R

######################################
### Load packages
######################################
library(tidyverse)
library(spatstat)

######################################
### Set parameters
######################################
which_sample = 1 #Specify which COCO sample to use (1, 2, ..., 77)
length_seq = 1000 #Length of estimated function
wd_data = "./data/coco_mw_samples/" #Working directory to MW samples 


######################################
### Load halo samples
######################################
cdm = read_csv(str_c(wd_data, "cdm_", which_sample, ".csv"))
wdm = read_csv(str_c(wd_data, "wdm_", which_sample, ".csv"))
xlim = range(cdm[,1], wdm[,1]) #Define boundary of sample region
ylim = range(cdm[,2], wdm[,2])
zlim = range(cdm[,3], wdm[,3])

######################################
### Compute G-Functions
######################################

### Define class as 3D point pattern
ppc = pp3(cdm$V1, cdm$V2, cdm$V3, 
           box3(xrange=xlim,yrange=ylim,zrange=zlim),
           marks=NULL)

ppw = pp3(wdm$V1, wdm$V2, wdm$V3, 
          box3(xrange=xlim,yrange=ylim,zrange=zlim),
          marks=NULL)
  
### G-function (CDM)
gc = G3est(ppc, correction="best", nrval = length_seq)
seq_gc = gc$r
fun_gc = gc$km

### G-function (WDM)
gw = G3est(ppw, correction="best", nrval = length_seq)
seq_gw = gw$r
fun_gw = gw$km


######################################
### Visualization of G-Functions
######################################
#Define data frame for plotting
df = tibble(tseq = c(seq_gc, seq_gw),
            gfun = c(fun_gc, fun_gw),
            type = c(rep("CDM",length_seq),
                     rep("WDM",length_seq)))

ggplot(df, aes(x=tseq, y=gfun, color=type, shape=type)) +
  geom_point() +
  geom_line() +
  xlab("Distance (Mpc)") +
  ylab(expression(F[G](t))) +
  labs(color=NULL, shape=NULL) +
  theme_bw() +
  theme(text = element_text(size=20))






