# DarkMatterTDA


This repository includes code, data, and summaries corresponding to the paper
"Uncovering Small-Scale Differences in the Large-Scale Structure of the Universe with Persistent Homology"

The code uses R (https://www.r-project.org/) and assumes you have installed the following packages:

sfsmisc 
spatstat
TDA
Tidyverse
ggpmisc (for plotting)
latex2exp (for plotting)
scales (for plotting)


Please reference using the following citation:
[CK2022] Cisewski-Kehe, J., Fasy, B.T., Hellwing, W., Lovell, M.R., Drozda, P. and Wu, M., 2022. Differentiating small-scale subhalo distributions in CDM and WDM models using persistent homology. arXiv preprint arXiv:2204.00443.


If you would like additional background on the math/statistics, please see the following:
Berry, E., Chen, Y.C., Cisewski-Kehe, J. and Fasy, B.T., 2020. Functional summaries of persistence diagrams. Journal of Applied and Computational Topology, 4(2), pp.211-262.


If you come across any issues or bugs, please contact Jessi Cisewski-Kehe at jjkehe@wisc.edu.


######################################## Data and Summaries

/data/coco_mw_samples
- The 77 CDM and WDM samples used in [CK2022]:  cdm_i.csv and wdm_i.csv for i=1,2, ..., 77
--77 samples of complete spatial randomness are also included as poisson_i.csv
- Each sample includes a MW-analog dark matter halo at the center, and other haloes within a 3 Mpc radius


/data/coco_mw_persistence_diagrams
- The Vietoris-Rips persistence diagrams computed for the CDM and WDM samples in /data/coco_mw_samples
--cdm_rips_i.txt and wdm_i.txt for i=1,2,...,77
- These were computed using the Ripser code by Bauer (2021) available at https://github.com/Ripser/ripser
- (Minor modifications were made to the Ripser code to produce output in the format of the included persistence diagrams)


/data/coco_mw_functional_summaries
- Includes the functional summaries computed for the persistence diagrams or the MW samples (in the case of the G-functions and 2PCF)
- The functions are saved as .csv files where column i is the functional summary corresponding to sample i
- The names of the files indicate the functional summary:
---landk = landscape functions for homology dimension k 
---ec = Euler characteristic functions
---eck = Betti functions for homology dimension k
---silk_p5 = silhouette functions for homology dimension k with p = 0.5
---silk_1 = silhouette functions for homology dimension k with p = 1
---silk_2 = silhouette functions for homology dimension k with p = 2
---pcf = two-point correlation functions (2PCF)
---gfun = G-functions
- The _cdm.csv = functions for the CDM samples; _wdm.csv = functions for the WDM samples
- The corresponding filtration parameter sequence for the persistence diagram-based functions use the following:
--- tda0_tseq.csv = sequence for H0
--- tda1_tseq.csv = sequence for H1
--- tda2_tseq.csv = sequence for H2
--- tda012_tseq.csv = sequence for the EC functions
- The corresponding filtration parameter sequence for the other functions are:
--- gfun_tseq.csv = sequence for G-functions
--- pcf_tseq.csv = sequence for 2PCF functions



######################################## Main Code

ComputeFunctionalSummaries.R
- Computes landscape functions, silhouette functions, and Euler characteristic and Betti functions from the persistence diagrams in /data/coco_mw_persistence_diagrams
- Each functional summary is saved individually as .rds files
- Note that you can also access pre-computed functional summaries in /data/coco_mw_functional_summaries where they are saved as .csv files in groups based on the functional summary and CDM/WDM status
- This code uses the getLandscape.R script found in folder r_functions


ComputePermutationPvalues.R
- This script contains code to compute traditional and matched-pairs permutation p-values between CDM and WDM using a specified functional summary as the test statistic
- Uses the pre-computed functional summaries in data/coco_mw_functional_summaries/
- User needs to select functional summary and tseq; homology dimension of functional summary and tseq need to match (see code for details)


######################################## Supplemental R functions
/r_functions

getLandscape.R
- Computes landscape functions
- Can also use `landscape()` from TDA package, but found the getLandscape.R functions to be faster


######################################## Additional Code
The previous code carries out the primary analysis of [CK2022]. Below is additional code for some graphics or supplementary analysis.


Confidence Band graphics (Figures 6 and 8 of [CK2022])
- R script:  ComputeConfidenceBands.R


MaxPersistence and Average Persistence Graphics (Figure 7 of [CK2022])
- R script:  MaxAvePersistenceGraphic.R


Test for complete spatial randomness (Figures 9 of [CK2022])
- R script:  CompleteSpatialRandomnessMCTest.R


Comparison method calculations
Calculation of the Persistence Diagram Test (PDT)
- R script:  ComputePDT.R


Calculation of the G-functions 
- R script:  ComputeGfunction.R  


Calculation of the two-point correlation function (2PCF)
- the Landy-Szalay estimator was used for these functions
- More description is to come (this was not computed in R)!
























