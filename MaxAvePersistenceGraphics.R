###This code produces the MaxPersistence and Average Persistence Graphics 

######################################
### Load packages
######################################
###The following R packages are necessary to run the code
library(tidyverse)
library(TDA)
library(scales) #muted()
library(ggpmisc) # puts text in top-right of plot
library(latex2exp) # Math notation in axis labels

######################################
### Set parameters and working directory
######################################
wd_data <- "./data/coco_mw_persistence_diagrams/" 
num_samp <- 77  #Number of regions sampled


######################################
### Define new functions
######################################
getPersistence <- function(Diag, tseq, which_dim){
  ##Returns the average persistence of homology group generators that persist
  ###at the values of tseq
  ##Diag = persistence diagram as (dimension, birth, death)
  ##tseq = grid to estimate average persistence
  ##which_dim = homology dimension of interest (0, 1, or 2)
  diag <- matrix(Diag[Diag[,1]==which_dim, 2:3], ncol = 2)
  out <- c()
  for(ii in 1:length(tseq)){
    get_features <- which(diag[,1]<=tseq[ii] & diag[,2]>= tseq[ii])
    if(length(get_features)>0){
      new_diag <- matrix(diag[get_features, ], ncol=2)
      out[ii] <- mean(new_diag[,2]-new_diag[,1])
    } else{out[ii] <- 0}
  }
  return(out)	
}

getMax <- function(Diag, tseq, which_dim){
  ##Returns the maximum persistence of homology group generators that persist
  ###at the values of tseq
  ##Diag = persistence diagram as (dimension, birth, death)
  ##tseq = grid to estimate average persistence
  ##which_dim = homology dimension of interest (0, 1, or 2)
  diag <- matrix(Diag[Diag[,1]==which_dim, 2:3], ncol = 2)
  out <- c()
  for(ii in 1:length(tseq)){
    get_features <- which(diag[,1]<=tseq[ii] & diag[,2]>= tseq[ii])
    if(length(get_features)>0){
      new_diag <- matrix(diag[get_features, ], ncol=2)
      out[ii] <- max(new_diag[,2]-new_diag[,1])
    } else{out[ii] <- 0}
  }
  return(out)	
}


######################################
### Load persistence diagrams
######################################
wdiag <- list()
cdiag <- list()
for(ii in 1:num_samp){
  w_temp <- read_delim(paste0(wd_data, "wdm_rips_",ii, ".txt"), 
                       delim=",", skip=1, col_names=FALSE, 
                       show_col_types=FALSE) %>%
    rename(dim=X1, birth=X2, death=X3) %>%
    arrange(dim, desc(death)) 
  wdiag[[ii]] <- as.matrix(w_temp, ncol=3)
  
  c_temp <- read_delim(paste0(wd_data, "cdm_rips_",ii, ".txt"), 
                       delim=",", skip=1, col_names=FALSE, 
                       show_col_types=FALSE) %>%
    rename(dim=X1, birth=X2, death=X3) %>%
    arrange(dim, desc(death)) 
  cdiag[[ii]] <- as.matrix(c_temp, ncol=3)
}


### Define tseq
tseq <- seq(0, 2.5, by=0.05)

######################################
### Compute summaries
######################################
### Compute summaries by homology group dimension
### Note: The H0 infinity point is removed
persistence_cdm0 <- sapply(1:num_samp, function(ii) getPersistence(cdiag[[ii]][-1,], tseq, 0))
persistence_wdm0 <- sapply(1:num_samp, function(ii) getPersistence(wdiag[[ii]][-1,], tseq, 0))
max_cdm0 <- sapply(1:num_samp, function(ii) getMax(cdiag[[ii]][-1,], tseq, 0))
max_wdm0 <- sapply(1:num_samp, function(ii) getMax(wdiag[[ii]][-1,], tseq, 0))
#
persistence_cdm1 <- sapply(1:num_samp, function(ii) getPersistence(cdiag[[ii]], tseq, 1))
persistence_wdm1 <- sapply(1:num_samp, function(ii) getPersistence(wdiag[[ii]], tseq, 1))
max_cdm1 <- sapply(1:num_samp, function(ii) getMax(cdiag[[ii]], tseq, 1))
max_wdm1 <- sapply(1:num_samp, function(ii) getMax(wdiag[[ii]], tseq, 1))
#
persistence_cdm2 <- sapply(1:num_samp, function(ii) getPersistence(cdiag[[ii]], tseq, 2))
persistence_wdm2 <- sapply(1:num_samp, function(ii) getPersistence(wdiag[[ii]], tseq, 2))
max_cdm2 <- sapply(1:num_samp, function(ii) getMax(cdiag[[ii]], tseq, 2))
max_wdm2 <- sapply(1:num_samp, function(ii) getMax(wdiag[[ii]], tseq, 2))


### Combine the functions into a single data frame for plotting:
df <- tibble(time = rep(tseq, 6),  
       dim = c(rep("0", 2*length(tseq)),
               rep("1", 2*length(tseq)),
               rep("2", 2*length(tseq))),
       fun = c(rep("AvePers", length(tseq)),
               rep("MaxPers", length(tseq)),
               rep("AvePers", length(tseq)),
               rep("MaxPers", length(tseq)),
               rep("AvePers", length(tseq)),
               rep("MaxPers", length(tseq))),
       y = c(apply(persistence_cdm0-persistence_wdm0, 1, mean),
             apply(max_cdm0-max_wdm0, 1, mean),
             apply(persistence_cdm1-persistence_wdm1, 1, mean),
             apply(max_cdm1-max_wdm1, 1, mean),
             apply(persistence_cdm2-persistence_wdm2, 1, mean),
             apply(max_cdm2-max_wdm2, 1, mean)),
       y_sd = c(apply(persistence_cdm0-persistence_wdm0, 1, sd)/sqrt(num_samp),
                apply(max_cdm0-max_wdm0, 1, sd)/sqrt(num_samp),
                apply(persistence_cdm1-persistence_wdm1, 1, sd)/sqrt(num_samp),
                apply(max_cdm1-max_wdm1, 1, sd)/sqrt(num_samp), 
                apply(persistence_cdm2-persistence_wdm2, 1, sd)/sqrt(num_samp),
                apply(max_cdm2-max_wdm2, 1, sd)/sqrt(num_samp))
)



###AvePersistence Graphic
ylab0 = TeX("$\\bar{F}_{diff, AvePers}$") #"AvePers"
df %>%
  filter(fun %in% c("AvePers")) %>%
  ggplot(aes(x=time, y=y, color=dim, shape=dim)) +
  geom_point(size=3) +
  geom_line(size=1) +
  geom_errorbar(aes(ymin=y-y_sd, ymax=y+y_sd), width=.01) +
  geom_vline(xintercept=seq(0, max(tseq), by=.1), 
             color="gray", linetype="dotted") +
  geom_hline(yintercept=0) +
  scale_color_discrete(name = NULL, 
                       labels = c(expression(H[0]), expression(H[1]), expression(H[2]))) +
  scale_shape_discrete(name = NULL, 
                       labels = c(expression(H[0]), expression(H[1]), expression(H[2]))) +
  theme_bw() +
  theme(text = element_text(size=40),
        legend.key.size = unit(.75, 'cm'),
        legend.title = element_blank(),
        legend.position=c(.85,.25),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black")) +
  labs(x = "Distance [Mpc]",
       y=ylab0) 



###MaxPersistence Graphic
ylab0 = TeX("$\\bar{F}_{diff, MaxPers}$") #"MaxPers"
df %>%
  filter(fun %in% c("MaxPers")) %>%
  ggplot(aes(x=time, y=y, color=dim, shape=dim)) +
  geom_point(size=3) +
  geom_line(size=1) +
  geom_errorbar(aes(ymin=y-y_sd, ymax=y+y_sd), width=.01) +
  geom_vline(xintercept=seq(0, max(tseq), by=.1), 
             color="gray", linetype="dotted") +
  geom_hline(yintercept=0) +
  scale_color_discrete(name = NULL, 
                       labels = c(expression(H[0]), expression(H[1]), expression(H[2]))) +
  scale_shape_discrete(name = NULL, 
                       labels = c(expression(H[0]), expression(H[1]), expression(H[2]))) +
  theme_bw() +
  theme(text = element_text(size=40),
        legend.key.size = unit(.75, 'cm'),
        legend.title = element_blank(),
        legend.position=c(.85,.25),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black")) +
  labs(x = "Distance [Mpc]",
       y=ylab0) 


