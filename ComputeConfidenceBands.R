###This code computes the produces global confidence bands for the functional
###summaries and plots them
###Note:  the paper uses 1000 bootstraps, but the code is set to use 100 to 
####reduce the computational burden

######################################
### Load packages
######################################
###The following R packages are necessary to run the code
library(tidyverse)
library(TDA)
library(scales) #muted()
library(ggpmisc) # puts text in top-right of plot
library(latex2exp) #for axis labels


######################################
### Set parameters and working directory
######################################
###The following sets up objects that are needed for computing the confidence
####bands for all the tests considered and for each homology dimension
###Note:  the Persistence Diagram Tests (PDT) test statistic is not a function
####and therefore does not have confidence bands computed
###Identifiers for the tests
tests = c("land0", "land1", "land2", 
           "sil0_p5", "sil0_1", "sil0_2", "sil1_p5", "sil1_1", "sil1_2", 
           "sil2_p5", "sil2_1", "sil2_2", 
           "ec", "ec0", "ec1", "ec2",
           "pdt0","pdt1","pdt2", #Note: PDT test statistic is not a function
           "gfun","pcf")


###The tseq for each test
test_tseq = c("tseq0", "tseq1", "tseq2", 
               "tseq0", "tseq0", "tseq0", 
               "tseq1", "tseq1", "tseq1", 
               "tseq2", "tseq2", "tseq2", 
               "tseq8", "tseq8", "tseq8", "tseq8",
               NA,NA,NA,NA,NA)


###The homology group dimension for each test
tests_hom = c(0,1,2,0,0,0,1,1,1,2,2,2,8,0,1,2,0,1,2,NA,NA)


###X axis labels for each test 
tests_x = c(rep("(Birth+Death)/2 (Mpc)", 12),
             rep("T (Mpc)",4),
             NA,NA,NA,
             rep("Distance (Mpc)",2))


###Y axis labels for each test 
tests_y = c(expression(F[land](1:10,0)),
             expression(F[land](1:10,1)),
             expression(F[land](1:10,2)),
             expression(F[sil](0,p==0.5)),
             expression(F[sil](0,p==1)),
             expression(F[sil](0,p==2)),
             expression(F[sil](1,p==0.5)),
             expression(F[sil](1,p==1)),
             expression(F[sil](1,p==2)),
             expression(F[sil](2,p==0.5)),
             expression(F[sil](2,p==1)),
             expression(F[sil](2,p==2)),
             expression(F[ec]),
             expression(F[betti](0)),
             expression(F[betti,cdm](1)-F[betti,wdm](1)),
             expression(F[betti](2)),
             NA, #PDT test statistic is not a function
             NA,
             NA,
             expression(F[G]),
             expression(F[PCF]))


###Y axis labels for each test (difference plots)
tests_y_diff = c(TeX("$\\bar{F}_{diff, land}(1-10,0)$"),
                  TeX("$\\bar{F}_{diff, land}(1-10,1)$"),
                  TeX("$\\bar{F}_{diff, land}(1-10,2)$"),
                  TeX("$\\bar{F}_{diff, sil}(0)$"),
                  TeX("$\\bar{F}_{diff, sil}(0)$"),
                  TeX("$\\bar{F}_{diff, sil}(0)$"),
                  TeX("$\\bar{F}_{diff, sil}(1)$"),
                  TeX("$\\bar{F}_{diff, sil}(1)$"),
                  TeX("$\\bar{F}_{diff, sil}(1)$"),
                  TeX("$\\bar{F}_{diff, sil}(2)$"),
                  TeX("$\\bar{F}_{diff, sil}(2)$"),
                  TeX("$\\bar{F}_{diff, sil}(2)$"),
                  TeX("$\\bar{F}_{diff, ec}$"),
                  TeX("$\\bar{F}_{diff, betti}(0)$"),
                  TeX("$\\bar{F}_{diff, betti}(1)$"),
                  TeX("$\\bar{F}_{diff, betti}(2)$"),
                  NA, #PDT test statistic is not a function
                  NA,
                  NA,
                  TeX("$\\bar{F}_{diff, G}$"),
                  TeX("$\\bar{F}_{diff, 2PCF}$"))


###Set working directory to functional summaries
wd_data = "./data/coco_mw_functional_summaries/" 


###Parameters
num_samp = 77  #Number of regions sampled
num_land = 10 #Number of landscape layers for test


###Load tseq by homology dimension
tseq0 = read.csv(paste0(wd_data,"tda0_tseq.csv"), header=FALSE)[[1]] #H0
tseq1 = read.csv(paste0(wd_data,"tda1_tseq.csv"), header=FALSE)[[1]] #H1
tseq2 = read.csv(paste0(wd_data,"tda2_tseq.csv"), header=FALSE)[[1]] #H2
tseq8 = read.csv(paste0(wd_data,"tda012_tseq.csv"), header=FALSE)[[1]] #EC/Betti functions
tseqg = read.csv(paste0(wd_data,"gfun_tseq.csv"), header=FALSE)[[1]] #G-function
tseq2pcf = read.csv(paste0(wd_data,"pcf_tseq.csv"), header=FALSE)[[1]] #2PCF


######################################
### Define new functions
######################################
#####----------------------------- Define Confidence Band
getConfBand = function(confidence_level=0.95, NBoot=100, test_functions){
  ### Computes a bootstrap confidence band for functions
  ###confidence_level = confidence level between 0 and 1 (defaults to 0.95)
  ###NBoot = number of bootstrap samples (defaults to 100)
  ###functions = each column is a functional summary
  nfunctions = ncol(test_functions)
  mean0 = apply(test_functions,1,mean)
  mean_functions = sapply(1:NBoot, function(ii) 
    abs(apply(test_functions[,sample(1:nfunctions, nfunctions, replace=TRUE)],1,mean)-mean0))
  max_functions = apply(mean_functions, 2, max)
  band_width = quantile(max_functions,confidence_level)
  
  df = data.frame(mean = mean0, lower=mean0-band_width, upper=mean0+band_width)
  return(df)
}

################################################################################
########################################Confidence Band Computations
################################################################################
###Get CDM and WDM confidence bands
df = tibble(tseq = NULL, type = NULL, test=NULL, dim=NULL,
             mean = NULL, lower = NULL,
             upper=NULL)

confidence_level = 0.95
Nboot = 100 #1000 is used in the paper
for(ii in 1:length(tests)){
  print(paste0("Test:  ", tests[ii], ", Number ", ii, " of ", length(tests)))
  which_test = tests[ii]
  test_index = which(tests==which_test)

  if(substr(which_test,1,3)!="pdt"){ #No PDT tests
  for(which_type in c("cdm","wdm")){ #Compute for CDM and WDM samples
    ###---------------------------landscapes functions
    if(substr(which_test,1,4)=="land"){
      tseq = get(test_tseq[test_index])
      fun = as.matrix(read.csv(str_c(wd_data,which_test,"_",which_type, ".csv"),
                            header=FALSE))
      tseq = seq(min(tseq),min(tseq)+num_land*length(tseq)*diff(tseq)[1],
                       length.out = num_land*length(tseq)) #maintain same spacing
    ###---------------------------G-functions
    }else if(which_test=="gfun"){
      tseq = tseqg
      fun = as.matrix(read.csv(str_c(wd_data,which_test,"_",which_type, ".csv"),
                                header=FALSE))
    ###---------------------------2PCF
      }else if(which_test=="pcf"){
        tseq = tseq2pcf
        fun = as.matrix(read.csv(str_c(wd_data,which_test,"_",which_type, ".csv"),
                                  header=FALSE))
    ###---------------------------Other functional summaries
        }else{#Other functions
          tseq = get(test_tseq[test_index])
          fun = as.matrix(read.csv(str_c(wd_data,which_test,"_",which_type, ".csv"),
                                    header=FALSE))
    }

    band_out = getConfBand(confidence_level, Nboot, fun)
    df0 = tibble(tseq = tseq, type = which_type, test=which_test, dim=tests_hom[test_index],
                  mean = band_out$mean, lower = band_out$lower,
                  upper=band_out$upper)
    df = bind_rows(df,df0)
  }}
}


###Bootstrap computations for the functional summary *differences*
df_diff = tibble(tseq = NULL, type = NULL, test=NULL, dim=NULL,
             mean = NULL, lower = NULL,
             upper=NULL)

confidence_level = 0.95
Nboot = 100 #1000 is used in the paper
for(ii in 1:length(tests)){
  print(paste0("Test:  ", tests[ii], ", Number ", ii, " of ", length(tests)))
  which_test = tests[ii]
  test_index = which(tests==which_test)
  
  if(substr(which_test,1,3)!="pdt"){ #No PDT tests
      ###---------------------------landscapes functions
      if(substr(which_test,1,4)=="land"){
        tseq = get(test_tseq[test_index])
        func = as.matrix(read.csv(str_c(wd_data,which_test,"_cdm.csv"),
                                  header=FALSE))
        funw = as.matrix(read.csv(str_c(wd_data,which_test,"_wdm.csv"),
                                   header=FALSE))
        tseq = seq(min(tseq),min(tseq)+num_land*length(tseq)*diff(tseq)[1],
                    length.out = num_land*length(tseq)) #maintain same spacing
        ###---------------------------G-functions
      }else if(which_test=="gfun"){
        tseq = tseqg
        func = as.matrix(read.csv(str_c(wd_data,which_test,"_cdm.csv"),
                                  header=FALSE))
        funw = as.matrix(read.csv(str_c(wd_data,which_test,"_wdm.csv"),
                                    header=FALSE))
        ###---------------------------2PCF
      }else if(which_test=="pcf"){
        tseq = tseq2pcf
        func = as.matrix(read.csv(str_c(wd_data,which_test,"_cdm.csv"),
                                  header=FALSE))
        funw = as.matrix(read.csv(str_c(wd_data,which_test,"_wdm.csv"),
                                   header=FALSE))
        ###---------------------------Other functional summaries
      }else{#Other functions
        tseq = get(test_tseq[test_index])
        func = as.matrix(read.csv(str_c(wd_data,which_test,"_cdm.csv"),
                                  header=FALSE))
        funw = as.matrix(read.csv(str_c(wd_data,which_test,"_wdm.csv"),
                                   header=FALSE))
      }
      fun = func - funw
      band_out = getConfBand(confidence_level, Nboot, fun)
      df0 = tibble(tseq = tseq, type = which_type, test=which_test, dim=tests_hom[test_index],
                    mean = band_out$mean, lower = band_out$lower,
                    upper=band_out$upper)
      df_diff = bind_rows(df_diff,df0)
    }
}

################################################################################
########################################Graphics
################################################################################
###Define color and line types:
colors = c("cdm" = "blue", "wdm" = "red")
lines = c("cdm" = "solid", "wdm" = "dotted")

###Confidence band plots for the functional summaries
####Specify the test index below
test_index = 15
(tests[test_index]) #Prints function type
df %>%
  filter(test==tests[test_index]) %>%
  ggplot(aes(x=tseq)) + 
  geom_line(aes(y=mean, color=type, linetype=type), size=2, show.legend=TRUE) + 
  geom_ribbon(aes(tseq,mean, ymin=lower, ymax=upper, color=type, fill=type), linetype=2, alpha=0.1) +
  labs(title = paste0("Test: ", tests[test_index]),
       x = tests_x[test_index],
       y = tests_y[test_index],
       color = NULL,
       linetype=NULL,
       fill=NULL) +
  scale_color_manual(values = colors,
                     labels = c("CDM","WDM"),
                     aesthetics = c("color", "fill")) +
  scale_linetype_manual(values = lines,
                        labels = c("CDM", "WDM")) +
  theme_bw() +
  theme(text = element_text(size=40),
        title = element_text(size=20),
        legend.title = element_blank(),
        legend.key.size = unit(.75, 'cm'),
        legend.position=c(.85,.9),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black")) 




###Confidence band plots for the functional summaries *differences*
####Specify the test index below
####See code below for (i) improved graphic for landscape functions, or
####(ii) including all silhouette parameters (p=0.5, 1, 2) on a single plot
test_index = 15
(tests[test_index]) #Prints function type
df_diff %>%
  filter(test==tests[test_index]) %>%
  ggplot(aes(x=tseq)) + 
  geom_line(aes(y=mean), size=2, show.legend=TRUE) + 
  geom_ribbon(aes(tseq,mean, ymin=lower, ymax=upper), linetype=2, alpha=0.1) +
  labs(title = NULL,
       x = tests_x[test_index],
       y = tests_y_diff[test_index]) +
  geom_hline(aes(yintercept=0)) +
  theme_bw() +
  theme(text = element_text(size=40),
        title = element_text(size=40),
        legend.title = element_blank(),
        legend.key.size = unit(1.75, 'cm'),
        legend.position=c(.75,.5),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black"))


###Confidence band plots for Landscape *differences*
####with an adjustment to the x-axis
test_index = 1 #Landscapes H0:1, H1:2, H2:3
(tests[test_index]) #Prints function type
tseq = df_diff %>% filter(test==tests[test_index]) %>% pull(tseq)
tseq_index = c(0,length(tseq)/num_land)
tseq_labels = c(0,tseq[length(tseq)/num_land])
df_diff %>%
  filter(test==tests[test_index]) %>%
  ggplot(aes(x=tseq)) + 
  geom_line(aes(y=mean), size=2, show.legend=TRUE) + 
  geom_ribbon(aes(tseq,mean, ymin=lower, ymax=upper), linetype=2, alpha=0.1) +
  labs(title = NULL,
       x = tests_x[test_index],
       y = tests_y_diff[test_index]) +
  geom_hline(aes(yintercept=0)) +
  theme_bw() +
  theme(text = element_text(size=40),
        title = element_text(size=40),
        legend.title = element_blank(),
        legend.key.size = unit(1.75, 'cm'),
        legend.position=c(.75,.5),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black"))+
  geom_vline(xintercept = c(0,tseq[c(1:10*1000)]), color="red", 
             linetype="dotted",size=2)+
  scale_x_continuous(breaks=c(min(tseq),tseq[tseq_index]), 
                   labels=round(tseq_labels,2))


###Confidence band plots for Silhouette *differences*
####plotted on same graphic
test_index = 4:6 #H0:4:6, H1:7:9, H2:10:12
(tests[test_index]) #Prints function type

if(min(test_index)==4){
  ylab0 = TeX("$\\bar{F}_{diff, sil}(0)$")
  colors = c("sil0_p5" = "red", "sil0_1" = "blue", "sil0_2"="green")
  lines = c("sil0_p5" = "solid", "sil0_1" = "dashed", "sil0_2"="dotted")
}else if(min(test_index)==7){
  ylab0 = TeX("$\\bar{F}_{diff, sil}(1)$")
  colors = c("sil1_p5" = "red", "sil1_1" = "blue", "sil1_2"="green")
  lines = c("sil1_p5" = "solid", "sil1_1" = "dashed", "sil1_2"="dotted")
}else if(min(test_index)==10){
  ylab0 = TeX("$\\bar{F}_{diff, sil}(2)$")
  colors = c("sil2_p5" = "red", "sil2_1" = "blue", "sil2_2"="green")
  lines = c("sil2_p5" = "solid", "sil2_1" = "dashed", "sil2_2"="dotted")
}

df_diff %>%
  filter(test%in%tests[test_index]) %>%
  ggplot(aes(x=tseq)) + 
  geom_line(aes(y=mean, color=test, linetype=test), size=2, show.legend=TRUE) + 
  geom_ribbon(aes(tseq,mean, ymin=lower, ymax=upper, color=test, fill=test), linetype=2, alpha=0.1) +
  labs(title = NULL, 
       x = tests_x[test_index[1]],
       y = ylab0) +
  geom_hline(aes(yintercept=0)) +
  scale_color_manual(values = colors,
                     labels = c("p=0.5","p=1","p=2"),
                     aesthetics = c("color", "fill")) +
  scale_linetype_manual(values = lines,
                        labels = c("p=0.5","p=1","p=2")) +
  theme_bw() +
  theme(text = element_text(size=40),
        title = element_text(size=40),
        legend.title = element_blank(),
        legend.key.size = unit(1.75, 'cm'),
        legend.position=c(.82,.25), 
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black")) 






