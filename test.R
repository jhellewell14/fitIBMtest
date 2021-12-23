suppressPackageStartupMessages(library(ggplot2))
library(malariasimulation)
library(malariaEquilibrium)
library(reshape2)

year <- 365
month <- 30

# Fake observed prevalence at 85 days after start
data <- 0.7

# MCMC stuff
parTab <- data.frame(values=c(365, 50, 500),
                     names=c("sim_length",
                             "starting_EIR",
                             "human_population"),
                     fixed=c(1,0,1),
                     lower_bound=c(-1000,0, 0),
                     upper_bound=c(1000,100, 1000),
                     steps=c(0.1,1, 1),
                     stringsAsFactors=FALSE)




## To test that it's working
posterior <- my_creation_function(parTab, data, my_prior)
print(posterior(parTab$values))

## Update the proposal step size every 1000 iterations (opt_freq) for the first 5000 iterations 
## (adaptive_period) to aim for an acceptance rate of 0.44 (popt). After the adaptive period, 
## run for a further 10000 (iterations) steps. Save every nth rows, where n is "thin" (ie. 1 here).
## Write to disk every 100 iterations (post thinning). Note that the adaptive period is also saved
## to disk
mcmcPars <- c("iterations"=100,"popt"=0.44,"opt_freq"=1000,
              "thin"=1,"adaptive_period"=25,"save_block"=50)

## The MCMC code uses the parameter table. Here, we should specify some random starting
## points in the "values" column.
startTab <- parTab
startTab$values[2] <- c(57)

output <- run_MCMC(parTab=startTab, data=data, mcmcPars=mcmcPars, filename="test", 
                   CREATE_POSTERIOR_FUNC=my_creation_function, mvrPars=NULL, 
                   PRIOR_FUNC = my_prior, OPT_TUNING=0.2)
