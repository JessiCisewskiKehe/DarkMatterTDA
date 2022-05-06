###This code computes a Monte Carlo global band under complete spatial
####randomness (CSR) to test if the CDM or WDM COCO samples are consistent
####with a 3D homogeneous spatial Poisson point process.
###See details in the paper for the interpretation of the result.


######################################
### Load packages
######################################
###The following R packages are necessary to run the code
library(tidyverse)
library(spatstat)


######################################
### Set parameters and working directory
######################################
which_sample = 1 #Specify which COCO sample to use (1, 2, ..., 77)
length_seq = 1000 #Length of estimated function
Nmc = 19 #Number of Monte Carlo samples for test
rad_samp = 3 #Radius (Mpc) of MW sample region
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
### Function to generate Poisson sample
######################################
pointSphere <- function(i, rad=3){
  ##i=iteration index
  ##rad=radius of sample region
  ##Return:  a halo randomly sampled in sphere of radius rad
  proposal <- runif(3, -rad, rad)
  dist0 <- sqrt(sum(proposal^2))
  while(dist0 > 3){
    proposal <- runif(3, -rad, rad)
    dist0 <- sqrt(sum(proposal^2))
  }
  return(proposal)
}


############################################
### Compute G-Functions for CDM/WDM sample
############################################

### Define class as 3D point pattern
ppc = pp3(cdm$V1, cdm$V2, cdm$V3, 
          box3(xrange=xlim,yrange=ylim,zrange=zlim),
          marks=NULL)

ppw = pp3(wdm$V1, wdm$V2, wdm$V3, 
          box3(xrange=xlim,yrange=ylim,zrange=zlim),
          marks=NULL)

### G-function (CDM)
gc = G3est(ppc, correction="best", nrval = length_seq)
fun_gc = gc$km
tseq = gc$r 

### G-function (WDM)
gw = G3est(ppw, correction="best", nrval = length_seq)
fun_gw = gw$km

### Note gc$r and gw$r should be equal
sum(abs(gc$r-gw$r)) == 0 #This should return TRUE


############################################
### Compute G-Functions under CSR
############################################
n = nrow(cdm) #Number of halos in sample
vol3 <- 4/3*pi*rad_samp^3 #Volume of sample region
lambda_hat <- n/vol3 #Estimated lambda (i.e., expected # halos per unit volume)

#Poisson G-function given lambda_hat
g_poisson_true <- 1-exp(-lambda_hat*4/3*pi*tseq^3)
  

############################################
### Compute global bound for MC Test
############################################
g_poisson_sample <- matrix(0, nrow=length_seq, ncol=Nmc)
lim_ppp <- rbind(rep(-rad_samp,3), rep(rad_samp,3))
for(jj in 1:Nmc){
  # Get Poisson sample
  ppp_sample <- t(sapply(1:n, function(ii) pointSphere(ii, rad_samp)))
  # Define class as 3D point pattern
  ppp <- pp3(ppp_sample[,1], ppp_sample[,2], ppp_sample[,3], 
             box3(xrange=lim_ppp[,1],yrange=lim_ppp[,2],zrange=lim_ppp[,3]),
             marks=NULL)
  # Estimate G-function for sample
  g_poisson_sample[,jj] <- G3est(ppp, correction="best", 
                                 nrval = length_seq, 
                                 rmax = max(tseq))$km
  }
poisson_difference <- g_poisson_sample - g_poisson_true
max_difference <- max(abs(poisson_difference))
  

######################################
### Visualization of G-Functions
######################################
#Define data frame for plotting
df <- tibble(tseq=rep(tseq,3),
             dm=c(rep("CDM", length_seq), 
                  rep("WDM", length_seq),
                  rep("CSR", length_seq)),
             lower=rep(g_poisson_true-max_difference,3),
             upper=rep(g_poisson_true+max_difference,3),
             obs=c(fun_gc, fun_gw, g_poisson_true))
  
colors <- c("CDM" = "blue", "WDM" = "red", "CSR"="black")
line_type <- c("CDM" = "dashed", "WDM" = "dotted", "CSR"="solid")

#Note the x-axis limit is adjusted in the graphic below
df %>%
  ggplot(aes(x=tseq, y=obs)) +
  geom_line(aes(y=obs, color=dm, linetype=dm), size=2) +
  geom_ribbon(aes(x=tseq, ymin=lower, ymax=upper),
              linetype=2, fill="gray", alpha=0.5, show.legend=FALSE) +
  theme_bw() +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = line_type) +
  theme(text=element_text(size=20)) +
  theme(legend.key.size = unit(2, 'cm'),
        legend.position=c(.75,.4)) +
  labs(x = "t (Mpc)",
       y = expression(F[G](t)),
       color=NULL,
       linetype=NULL) +
  xlim(0,1.25) 

  
  
  