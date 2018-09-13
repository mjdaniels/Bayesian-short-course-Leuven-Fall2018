#only for the first run

install.packages("devtools")

install.packages(c("MCMCpack", "tidyverse", "gridExtra", "Rcpp", "RcppArmadillo", "mitools", "bindata","rjags","norm","runjags"))

library(devtools)
install_github("theodds/NiNBayes")


library(bindata)
library(norm)
library(runjags)
library(NiNBayes)

# code for schizo data

#dataset
schizo_data <- readRDS("schizo.rds")

#treatment arm - scaling the data (for the priors)
# here, location/scale transformation 
Y <- scale(schizo_data[[1]])

#control arm
#Y <- scale(schizo_data[[3]])

#create missing data indicator
R <- ifelse(is.na(Y), 0, 1)

# adjust to match data example slides - nrow (value 100 in the below) is the number of posterior samples kept and then used for G-computation
# since response is scaled, corresponds to a a prior which is uniform (0,sigma_j)
sens.param <- matrix(runif(100 * ncol(Y)), nrow = 100)

# obtain the posterior sample of the parameters in the working model
# K is the number of clusters/components
# n.save is the number of posterior samples to use for G-computation 
foo <- FitDPMN(Y,R, K = 5, n.adapt = 100, n.save = 100, thin = 1, n.chains = 1)

# G-computation step
# defaults to 10K iterations for each posterior sample 
bar <- GCompDPMN(foo$mcmc, Y, K = 5, sens.param = sens.param)

# output from GCompDPMN
#> names(bar)
#[1] "means"          "obs.means"      "drop.probs"     "haz"           
#[5] "cum.drop.probs"

# we want 'means' - these are the 'standardized' means at each time
# other outputs are for model checking - see the NINBayes package for details 

#rescale means at time J for treatment group
rescale.mean.timeJ = bar$means[,6]*sd(schizo_data[[1]][,6],na.rm=TRUE)+ mean(schizo_data[[1]][,6],na.rm=TRUE) 
# posterior mean of the mean outcome at time J
print(mean(rescale.mean.timeJ))
# 95% CI for the mean outcome at time J
print(quantile(rescale.mean.timeJ,c(.025,.975)))


# there are more complex possibilities in the NINBayes repository including the non-monotone analyses 