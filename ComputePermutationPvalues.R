######################################
### Overview
######################################
## This script contains code to compute traditional and matched-pairs permutation
###...p-values between CDM and WDM using a specified functional summary as the
###...test statistic
## Uses the pre-computed functional summaries in data/coco_mw_functional_summaries/
## User needs to select functional summary and tseq below under "Select function and tseq"

######################################
### Load packages
######################################
library(tidyverse)
library(TDA)
library(sfsmisc) #integrate.xy


######################################
### Set parameters and working directory
######################################
num_land <- 10 #Number of landscape layers to retain
num_samp <- 77 #Number of regions sampled
nPerm <- 100 #Number of permutations.  In the paper, 20000 were used.
wd_data <- "data/coco_mw_functional_summaries/" #Working directory to data files (the functional summaries)

######################################
### Select function and tseq
######################################
### Lists the function files and the tseq options
function_files <- list.files(wd_data)

### CAUTION:  make sure the tseq homology dimension matches the functional summary
####...homology dimension.  The example below uses 1-Betti functions
### Options for tseq:  
#### H0:  tda0_tseq.csv
#### H1:  tda1_tseq.csv
#### H2:  tda2_tseq.csv
#### EC function: tda012_tseq.csv
#### G-function:  gfun_tseq.csv
#### 2PCF:  pcf_tseq.csv
tseq_file_name = "tda1_tseq.csv"
tseq <- read.csv(paste0(wd_data,tseq_file_name), header=FALSE)$V1

###Select functional summary, make sure it matches the homology dimension of tseq
which_test <- "ec1"
fun_cdm <- read.csv(paste0(wd_data,which_test,"_cdm.csv"), header=FALSE)
fun_wdm <- read.csv(paste0(wd_data,which_test,"_wdm.csv"), header=FALSE)

###Landscapes:  If functional summaries are landscapes, need to adjust tseq to account for 
####...concatenation of landscape layers
if(substr(which_test,1,4)=="land"){#landscapes functions
  tseq <- seq(min(tseq),min(tseq)+num_land*length(tseq)*diff(tseq)[1],
              length.out = num_land*length(tseq)) #maintain same spacing
  }


######################################
### Define functions
######################################
funPerm <- function(i, fun_use, tseq){
  #Permutation test for functions
  #i = arbitrary, used for iteration
  #fun_use = matrix of both groups of functions; each column is a function
  ##First 77 columns are from group 1, next 77 columns are from group 2
  ##(code below assumes equal number of functions from each sample)
  #tseq = x-axis for functional summaries
  which_fun <- sample(1:ncol(fun_use), ncol(fun_use)/2, replace = FALSE)
  mean0 <- apply(fun_use[,which_fun],1,mean)		
  mean1 <- apply(fun_use[,-which_fun],1,mean)
  return(integrate.xy(tseq,abs(mean0 - mean1)))
}

funPermMatched <- function(i, fun_use1, fun_use2, tseq){
  # A Matched-Pairs Permutation test for functions
  #i = arbitrary, used for iteration
  #fun_use1 = matrix of functions from group 1; each column is a function
  #fun_use2 = matrix of functions from group 2; each column is a function
  ##(code below assumes functions with matching column indexes from fun_use1 and fun_use2
  ##are the matched samples...column i from fun_use1 is matched with column i from fun_use2)
  #tseq = x-axis for functional summaries
  ncols <- ncol(fun_use1)
  groups <- sample(1:2, ncols, replace = TRUE)
  which1 <- which(groups == 1)
  which2 <- which(groups == 2)
  
  fun_use_11 <- cbind(fun_use1[, which1], fun_use2[, which2])
  fun_use_22 <- cbind(fun_use1[, which2], fun_use2[, which1])
  mean0 <- apply(fun_use_11,1,mean)		
  mean1 <- apply(fun_use_22,1,mean)
  return(integrate.xy(tseq,abs(mean0 - mean1)))
}



######################################
### Compute p-values
######################################
print(paste0("Permutation Tests"))
###---------------------------------------------------------------Perform test

#Calculate mean functional summary for each group		
mean_cdm <- apply(fun_cdm,1,mean)		
mean_wdm <- apply(fun_wdm,1,mean)
  
#Calculate the test statistic
test_stat <- integrate.xy(tseq,abs(mean_cdm - mean_wdm))
  
#Combine the functions; run the permutation tests
fun_use <- cbind(fun_cdm, fun_wdm)
permDist <- unlist(lapply(rep(1,nPerm), funPerm, fun_use, tseq))
pval <- length(which(test_stat <= permDist))/nPerm
print(paste0("Traditional permutation p-value: ", pval, " with ", nPerm, " permutations"))

#Run the matched permutation tests
fun_use1 <- fun_cdm
fun_use2 <- fun_wdm
permDist_matched <- unlist(lapply(rep(1,nPerm), funPermMatched, fun_use1, fun_use2, tseq))
pval_matched <- length(which(test_stat <= permDist_matched))/nPerm
print(paste0("Matched-pairs permutation p-value: ", pval_matched, " with ", nPerm, " permutations"))











