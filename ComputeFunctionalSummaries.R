######################################
### Overview
######################################
## This script contains code to compute the functional summaries of persistence diagrams
## These functions are already available in /data/coco_mw_functional_summaries for the
###...persistence diagrams in /data/coco_mw_persistence_diagrams
## Each function type is save for individual samples

######################################
### Load packages
######################################
library(tidyverse)
library(TDA)

######################################
### Load functions
######################################
source("r_functions/getLandscape.R")

######################################
### Set parameters and working directory
######################################
num_land <- 10 #Number of landscape layers to retain
num_samp <- 77 #Number of regions sampled
wd_data <- "data/" #Working directory to data files
wd_save <- "data/" #Working directory to save functional summaries

######################################
### Load persistence diagrams
######################################
diagw <- list() # To store the WDM diagrams
diagc <- list() # To store the CDM diagrams
for(ii in 1:num_samp){
  print(paste("Diagram ", ii, " of ", num_samp))
  diagw[[ii]] <- read_delim(paste0(wd_data, "coco_mw_persistence_diagrams/wdm_rips_",ii, ".txt"), 
                                   delim=",", skip=1, col_names=FALSE,
                            show_col_types = FALSE) %>%
    rename(dim=X1, birth=X2, death=X3) %>%
    arrange(dim, desc(death))
  diagc[[ii]] <- read_delim(paste0(wd_data, "coco_mw_persistence_diagrams/cdm_rips_",ii, ".txt"), 
                            delim=",", skip=1, col_names=FALSE,
                            show_col_types = FALSE) %>%
    rename(dim=X1, birth=X2, death=X3) %>%
    arrange(dim, desc(death))
}


##############################################
### Define the filtration parameter sequence
##############################################
### These can also be loaded directly from /data/coco_mw_functional_summaries
### The code below find the range of all the H0 features, all the H1 features, etc.
####...in order to define the tseq (the horizontal axis for the functional summaries)
####...so each type of functional summary is comparable across all the diagrams
#### Note that the H0 infinity points are removed because they are assigned arbitrary values

####Get filtration parameter ranges for the CDM and WDM persistence diagrams by homology group dimension
diagc0 <- lapply(1:num_samp,function(ii) diagc[[ii]][diagc[[ii]][,1]==0,])
diagc0_range <- sapply(1:num_samp,function(ii) range(diagc0[[ii]][-1,2:3])) #Removes infinity point
diagc1_range <- sapply(1:num_samp,function(ii) range(diagc[[ii]][diagc[[ii]][,1]==1,2:3]))
diagc2_range <- sapply(1:num_samp,function(ii) range(diagc[[ii]][diagc[[ii]][,1]==2,2:3]))

diagw0 <- lapply(1:num_samp,function(ii) diagw[[ii]][diagw[[ii]][,1]==0,])
diagw0_range <- sapply(1:num_samp,function(ii) range(diagw0[[ii]][-1,2:3])) #Removes infinity point
diagw1_range <- sapply(1:num_samp,function(ii) range(diagw[[ii]][diagw[[ii]][,1]==1,2:3]))
diagw2_range <- sapply(1:num_samp,function(ii) range(diagw[[ii]][diagw[[ii]][,1]==2,2:3]))

####Get range by homology group dimension
range0 <- range(diagw0_range,diagc0_range)
range1 <- range(diagw1_range,diagc1_range)
range2 <- range(diagw2_range,diagc2_range) 
range012 <- range(range0,range1,range2)

####Define the functional summary grids based on ranges by homology group dimension
tseq0 <- seq(range0[1], range0[2], length.out = 1000) #grid for H0 summaries
tseq1 <- seq(range1[1], range1[2], length.out = 1000) #grid for H1 summaries
tseq2 <- seq(range2[1], range2[2], length.out = 1000) #grid for H2 summaries
tseq012 <- seq(range012[1], range012[2], length.out = 1000) #grid for EC functions 



##############################################
### Estimate the function summaries
##############################################
### The code below goes through each sample and computes functional summaries for
####...landscape functions, silhouettes (w/p=0.5, 1, 2), Betti functions, and
####...the Euler characteristic function
### These functions are saved in the directory `wd_save` as .rds files
for(i in 1:num_samp){#estimate the functional summaries
  # print(i)
  
  print(paste0("Landscapes for sample: ", i))
  #-----------------------------Landscape functions
  
  ###CDM birth and death times as matrix
  subdiagc0 <- as.matrix(diagc[[i]][diagc[[i]][,1]==0,2:3], ncol=2)#H0 birth, death
  subdiagc1 <- as.matrix(diagc[[i]][diagc[[i]][,1]==1,2:3],ncol=2)#H1 birth, death
  subdiagc2 <- as.matrix(diagc[[i]][diagc[[i]][,1]==2,2:3],ncol=2)#H2 birth, death
  
  ###H0 landscape for CDM diagrams
  cland0 <- getLandscape(subdiagc0[-1,], tseq0) #remove point at infinity
  
  ###H1 landscape for CDM diagrams
  if(length(subdiagc1)>0){ ###Check that H1 features exist, if not assign 0's
    cland1 <- getLandscape(subdiagc1, tseq1) 
  }else{cland1 <- matrix(0, nrow = length(tseq1), ncol = 1)}
  
  if(length(subdiagc2)>0){ ###Check that H2 features exist, if not assign 0's
    cland2 <- matrix(getLandscape(subdiagc2, tseq2), nrow = length(tseq2))
  }else{cland2 <- matrix(0, nrow = length(tseq2), ncol = 1)}
  
  saveRDS(cland0[,1:min(num_land,ncol(cland0))], paste0(wd_save,"land0_cdm_", i, ".rds"))
  saveRDS(cland1[,1:min(num_land,ncol(cland1))], paste0(wd_save,"land1_cdm_", i, ".rds"))
  saveRDS(cland2[,1:min(num_land,ncol(cland2))], paste0(wd_save,"land2_cdm_", i, ".rds"))
  
  
  ###Repeat above for WDM diagrams
  subdiagw0 <- as.matrix(diagw[[i]][diagw[[i]][,1]==0,2:3],ncol=2)
  subdiagw1 <- as.matrix(diagw[[i]][diagw[[i]][,1]==1,2:3],ncol=2)
  subdiagw2 <- as.matrix(diagw[[i]][diagw[[i]][,1]==2,2:3],ncol=2)
  wland0 = getLandscape(subdiagw0[-1,], tseq0) #remove point at infinity
  
  if(length(subdiagw1)>0){
    wland1 <- getLandscape(subdiagw1, tseq1) 
  }else{wland1 <- matrix(0, nrow = length(tseq1), ncol = 1)}
  
  if(length(subdiagw2)>0){
    wland2 <- getLandscape(subdiagw2, tseq2) 
  }else{wland2 <- matrix(0, nrow = length(tseq2), ncol = 1)}
  
  saveRDS(wland0[,1:min(num_land,ncol(wland0))], 
          paste0(wd_save,"land0_wdm_", i, ".rds"))
  saveRDS(wland1[,1:min(num_land,ncol(wland1))], 
          paste0(wd_save,"land1_wdm_", i, ".rds"))
  saveRDS(wland2[,1:min(num_land,ncol(wland2))], 
          paste0(wd_save,"land2_wdm_", i, ".rds"))
  

  print(paste0("Silhouettes for sample: ", i))
  #-----------------------------Silhouette functions (p = 0.5, 1, 2)
  #Remove the point at infinity for H0
  csil0_p5 <- silhouette(as.matrix(diagc[[i]][-1,],ncol=3), p = .5,  dimension = 0, tseq0)
  csil0_1 <- silhouette(as.matrix(diagc[[i]][-1,],ncol=3), p = 1,  dimension = 0, tseq0)
  csil0_2 <- silhouette(as.matrix(diagc[[i]][-1,],ncol=3), p = 2,  dimension = 0, tseq0)
  csil1_p5 <- silhouette(as.matrix(diagc[[i]],ncol=3), p = .5,  dimension = 1, tseq1)
  csil1_1 <- silhouette(as.matrix(diagc[[i]],ncol=3), p = 1,  dimension = 1, tseq1)
  csil1_2 <- silhouette(as.matrix(diagc[[i]],ncol=3), p = 2,  dimension = 1, tseq1)
  csil2_p5 <- silhouette(as.matrix(diagc[[i]],ncol=3), p = .5,  dimension = 2, tseq2)
  csil2_1 <- silhouette(as.matrix(diagc[[i]],ncol=3), p = 1,  dimension = 2, tseq2)
  csil2_2 <- silhouette(as.matrix(diagc[[i]],ncol=3), p = 2,  dimension = 2, tseq2)
  saveRDS(csil0_p5, paste0(wd_save,"sil0_p5_cdm_",i,".rds"))
  saveRDS(csil0_1, paste0(wd_save,"sil0_1_cdm_",i,".rds"))
  saveRDS(csil0_2, paste0(wd_save,"sil0_2_cdm_",i,".rds"))
  saveRDS(csil1_p5, paste0(wd_save,"sil1_p5_cdm_",i,".rds"))
  saveRDS(csil1_1, paste0(wd_save,"sil1_1_cdm_",i,".rds"))
  saveRDS(csil1_2, paste0(wd_save,"sil1_2_cdm_",i,".rds"))
  saveRDS(csil2_p5, paste0(wd_save,"sil2_p5_cdm_",i,".rds"))
  saveRDS(csil2_1, paste0(wd_save,"sil2_1_cdm_",i,".rds"))
  saveRDS(csil2_2, paste0(wd_save,"sil2_2_cdm_",i,".rds"))
  
  wsil0_p5 <- silhouette(as.matrix(diagw[[i]][-1,],ncol=3), p = .5,  dimension = 0, tseq0)
  wsil0_1 <- silhouette(as.matrix(diagw[[i]][-1,],ncol=3), p = 1,  dimension = 0, tseq0)
  wsil0_2 <- silhouette(as.matrix(diagw[[i]][-1,],ncol=3), p = 2,  dimension = 0, tseq0)
  wsil1_p5 <- silhouette(as.matrix(diagw[[i]],ncol=3), p = .5,  dimension = 1, tseq1)
  wsil1_1 <- silhouette(as.matrix(diagw[[i]],ncol=3), p = 1,  dimension = 1, tseq1)
  wsil1_2 <- silhouette(as.matrix(diagw[[i]],ncol=3), p = 2,  dimension = 1, tseq1)
  wsil2_p5 <- silhouette(as.matrix(diagw[[i]],ncol=3), p = .5,  dimension = 2, tseq2)
  wsil2_1 <- silhouette(as.matrix(diagw[[i]],ncol=3), p = 1,  dimension = 2, tseq2)
  wsil2_2 <- silhouette(as.matrix(diagw[[i]],ncol=3), p = 2,  dimension = 2, tseq2)
  saveRDS(wsil0_p5, paste0(wd_save,"sil0_p5_wdm_",i,".rds"))
  saveRDS(wsil0_1, paste0(wd_save,"sil0_1_wdm_",i,".rds"))
  saveRDS(wsil0_2, paste0(wd_save,"sil0_2_wdm_",i,".rds"))
  saveRDS(wsil1_p5, paste0(wd_save,"sil1_p5_wdm_",i,".rds"))
  saveRDS(wsil1_1, paste0(wd_save,"sil1_1_wdm_",i,".rds"))
  saveRDS(wsil1_2, paste0(wd_save,"sil1_2_wdm_",i,".rds"))
  saveRDS(wsil2_p5, paste0(wd_save,"sil2_p5_wdm_",i,".rds"))
  saveRDS(wsil2_1, paste0(wd_save,"sil2_1_wdm_",i,".rds"))
  saveRDS(wsil2_2, paste0(wd_save,"sil2_2_wdm_",i,".rds"))
  

  print(paste0("EC and Bettis for sample: ", i))
  #-----------------------------Euler characteristic and Betti functions
  #First column of "ec" is the EC, second column = Betti-0, third column = Betti - 1, etc.
  Diag <- as.matrix(diagc[[i]], ncol=3)
  maxdimension <- 2
  tseq <- tseq012
  betti <- matrix(NA, nrow = length(tseq), ncol = maxdimension + 1)
  ec <- matrix(0, nrow = length(tseq), ncol = 1)
  for(jj in 0:maxdimension){ 
    diag <- matrix(Diag[Diag[,1]==jj, 2:3], ncol = 2)
    betti[,jj+1] <- sapply(tseq, function(ii) length(which(diag[,1]<=ii & diag[,2]>= ii)))
    ec <- ec + (-1)^jj*(betti[,jj+1])
  }
  cec <- cbind(ec, betti)
  
  saveRDS(cec[,1], paste0(wd_save,"ec_cdm_",i,".rds"))
  saveRDS(cec[,2], paste0(wd_save,"ec0_cdm_",i,".rds"))
  saveRDS(cec[,3], paste0(wd_save,"ec1_cdm_",i,".rds"))
  saveRDS(cec[,4], paste0(wd_save,"ec2_cdm_",i,".rds"))
  
  Diag <- as.matrix(diagw[[i]], ncol=3)
  betti <- matrix(NA, nrow = length(tseq), ncol = maxdimension + 1)
  ec <- matrix(0, nrow = length(tseq), ncol = 1)
  for(jj in 0:maxdimension){ 
    diag <- matrix(Diag[Diag[,1]==jj, 2:3], ncol = 2)
    betti[,jj+1] <- sapply(tseq, function(ii) length(which(diag[,1]<=ii & diag[,2]>= ii)))
    ec <- ec + (-1)^jj*(betti[,jj+1])
  }
  wec <- cbind(ec, betti)
  saveRDS(wec[,1], paste0(wd_save,"ec_wdm_",i,".rds"))
  saveRDS(wec[,2], paste0(wd_save,"ec0_wdm_",i,".rds"))
  saveRDS(wec[,3], paste0(wd_save,"ec1_wdm_",i,".rds"))
  saveRDS(wec[,4], paste0(wd_save,"ec2_wdm_",i,".rds"))
}
