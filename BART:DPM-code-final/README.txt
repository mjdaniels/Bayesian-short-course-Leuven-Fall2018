BNP approach for causal inference on quantiles

The first seven lines of the R code, EHR_example-course.r show the packages that need to be installed including package from my GitHub site, Sequential BART

There are two steps to implement the proposed approach (BART_DPM).
1. The first step estimates propensity score using BART. 
   If no missing data, use R package 'BayesTree'. 
   If covariates have missing data, use Sequential BART approach (https://github.com/mjdaniels/SequentialBART).
   WE HAVE MISSING COVARIATES SO WE WILL USE SEQUENTIAL BART

2. The second step estimates marginal distribution of potential outcomes using DPM.
   Use modified function 'DPcdensity' in DPpakage.  However, this package needs to be modified to use.  The following steps can be done to do this: 

(a) download the package source (DPpackage version 1.1.6) from 
   https://cran.r-project.org/src/contrib/Archive/DPpackage/
   Unzip it and replace the DPcdensity.R file in the R folder and the DPdensityreg.f file in the src folder with the modified code found in BART:DPM-code-final/quantile_codes on my GitHub page

                                                   
(b) replace DPcdensity.r in the R folder and DPdensityreg.f in the src folder;
                                                   
(c) install the revised DPpackage. If Rstudio is used, open file->new project->Existing Directory, find the directory where DPpackage is saved (.../DPpackage_1.1-6/DPpackage), then click 'create project'. And build the package from Build-> Clean and Rebuild. The package will be installed. 


The output save.state$denspm (a matrix of nsave*ngrid) has pdf of potential outcomes and save.state$denspl has cdf of potential outcomes at the grid points.

We will illustrate with the EHR data (with noise) discussed during the course (EHR_data.cvs)
