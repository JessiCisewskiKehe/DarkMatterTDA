###This code estimates the Persistence Diagram Test (PDT) p-values
###Note: computing the bottleneck distances between the diagrams is computationally
####intensive.  You may consider running those computations on a cluster.

###This implementation is based off of the test proposed in
######Robinson, A. and Turner, K., 2017. Hypothesis testing for topological data analysis. 
######Journal of Applied and Computational Topology, 1(2), pp.241-261.

######################################
### Load packages
######################################
library(TDA)
library(tidyverse)

######################################
### Set parameters
######################################
which_hom = 1 #Specify the homology dimension
qq = 1 #Parameter in PDT test statistic
num_samp <- 77 #Number of regions sampled
nPerm <- 100 #Number of permutations.  In the paper, 20000 were used.
wd_data <- "./data/coco_mw_persistence_diagrams/" #Working directory to persistence diagrams 



######################################
### Define new functions
######################################
funPerm <- function(i, p0, p1, qq=1, DIST){
  #Permutation test for PDT
  #i = arbitrary, used for iteration
  #p0 = number of observations in sample 0
  #p1 = number of observations in sample 1
  #qq = exoponent for test statistic computation (defaults to 1)
  #DIST = distance matrix of bottleneck distances between samples 0 and 1
	which.fun0 <- sort(sample(1:(p0+p1), p0, replace = FALSE)) #Random sample for group 0
	which.fun1 <- 1:(p0+p1) #All sample indexes
	which.fun1 <- which.fun1[-which.fun0] #Remove random group 0 indexes
	
	#Compute test statistic
	test.stat1 <- sum((DIST[which.fun0, which.fun0])^qq)/(p0*(p0-1))+
	  sum((DIST[which.fun1, which.fun1])^qq)/(p1*(p1-1))
	return(test.stat1)
}

funPermMatched <- function(i, p0, p1, qq, DIST){
  # A Matched-Pairs Permutation test for PDT
  #i = arbitrary, used for iteration
  #p0 = number of observations in sample 0
  #p1 = number of observations in sample 1
  #qq = exoponent for test statistic computation (defaults to 1)
  #DIST = distance matrix of bottleneck distances between samples 0 and 
  groups <- sample(0:1, p0, replace = TRUE)
  # 0: go into group 0
  # 1: go into group 1
  which0 <- which(groups == 0) #Index for Group 0
  which1 <- which(groups == 1) #Index for Group 1
  
  which.fun0 <- c(which0, which1+num_samp)
  which.fun1 <- c(which1, which0+num_samp)
  
  test.stat1 <- sum((DIST[which.fun0, which.fun0])^qq)/(p0*(p0-1))+
    sum((DIST[which.fun1, which.fun1])^qq)/(p1*(p1-1))
  return(test.stat1)
}


######################################
### Load persistence diagrams
######################################
diagw <- list() # To store the WDM diagrams
diagc <- list() # To store the CDM diagrams
for(ii in 1:num_samp){
  print(paste("Diagram ", ii, " of ", num_samp))
  diagw[[ii]] <- read_delim(paste0(wd_data, "wdm_rips_",ii, ".txt"), 
                            delim=",", skip=1, col_names=FALSE,
                            show_col_types = FALSE) %>%
    rename(dim=X1, birth=X2, death=X3) %>%
    arrange(dim, desc(death))
  diagc[[ii]] <- read_delim(paste0(wd_data, "cdm_rips_",ii, ".txt"), 
                            delim=",", skip=1, col_names=FALSE,
                            show_col_types = FALSE) %>%
    rename(dim=X1, birth=X2, death=X3) %>%
    arrange(dim, desc(death))
}
pdiag_01 <- c(diagw, diagc)


######################################
###Permutation test computations
######################################
distance0 <- matrix(0, ncol = 2*num_samp, nrow = 2*num_samp)

###Compute the Bottleneck distances between all pairs of persistence diagrams
###Note:  these computations are slow, but can easily be run in parallel or on
####a computing cluster.
print(paste0("Homology group: ", which_hom))
for(ii in 1:(num_samp*2-1)){
  print(ii)
  distance0[(ii+1):length(pdiag_01),ii] <- sapply((ii+1):length(pdiag_01), 
                                                  function(jj) bottleneck(as.matrix(pdiag_01[[ii]],ncol=3), 
                                                                          as.matrix(pdiag_01[[jj]],ncol=3), 
                                                                          dimension = which_hom))
  }

### Calculate observed statistic
test_stat <- sum((distance0[1:num_samp, 1:num_samp])^qq)/(num_samp*(num_samp-1)) +
    sum((distance0[(num_samp+1):(2*num_samp), 
                   (num_samp+1):(2*num_samp)])^qq)/(num_samp*(num_samp-1))
print(paste0("Observed test statistic: ", test_stat))

### Run permutation test	
permDist <- sapply(rep(1,nPerm), function(ii) funPerm(ii, 
                                                      p0=num_samp,
                                                      p1=num_samp, 
                                                      qq, 
                                                      DIST = distance0))
  
### Output:  p-value and test statistic
pval <- length(which(test_stat >= permDist))/nPerm
print(paste0("Permutation p-value: ", pval))
  
### Run matched permutation test	
permDist_matched <- sapply(rep(1,nPerm), function(ii) funPermMatched(ii, 
                                                                     p0=num_samp, 
                                                                     p1=num_samp,
                                                                     qq,
                                                                     DIST = distance0))
  
### Output matched-pairs:  p-value (test statistic is same as above)
pval_matched <- length(which(test_stat >= permDist_matched))/nPerm
print(paste0("Permutation p-value (matched-pairs): ", pval_matched))








