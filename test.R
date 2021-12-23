suppressPackageStartupMessages(library(ggplot2))
library(malariasimulation)
library(malariaEquilibrium)
library(lazymcmc)
source("aux_funcs.R")

year <- 365
month <- 30
human_population <- 1000
starting_EIR <- 50
sim_length <- 90

## "Real" run to recreate
r_simparams <- get_parameters(
  list(
    human_population = human_population,
    model_seasonality = TRUE, # Let's try a bi-modal model
    g0 = 0.28605,
    g = c(0.20636, -0.0740318, -0.0009293),
    h = c(0.173743, -0.0730962, -0.116019),
    severe_prevalence_rendering_min_ages = 2*year,
    severe_prevalence_rendering_max_ages = 10*year,
    severe_enabled = 1
  )
)

r_simparams <- set_equilibrium(r_simparams, starting_EIR)

r_output <- run_simulation(sim_length, r_simparams)
r_prev <- (r_output$n_detect_730_3650 / r_output$n_730_3650)[85]

## Set up MCMC

# Fake observed prevalence at 85 days after start
data <- r_prev

# MCMC stuff
parTab <- data.frame(values=c(sim_length, starting_EIR, human_population),
                     names=c("sim_length",
                             "starting_EIR",
                             "human_population"),
                     fixed=c(1,0,1),
                     lower_bound=c(-1000,0, 0),
                     upper_bound=c(1000,100, 1000),
                     steps=c(0.1, 1, 1),
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
startTab$values[2] <- rnorm(1, starting_EIR, 5)

output <- run_MCMC(parTab=startTab, data=data, mcmcPars=mcmcPars, filename="test", 
                   CREATE_POSTERIOR_FUNC=my_creation_function, mvrPars=NULL, 
                   PRIOR_FUNC = my_prior, OPT_TUNING=0.2)

chain <- read.csv(output$file)
plot(coda::as.mcmc(chain[chain$sampno > mcmcPars["adaptive_period"],]))
