
## You MUST specify the arguments to this function as parTab, data then PRIOR_FUNC. 
## Use the `...` to pass additional arguments through.
my_creation_function <- function(parTab, data, PRIOR_FUNC, ...){
  ##############################
  ## This is where you would manipulate all
  ## of your model arguments. For example,
  ## unpackaging the stuff from `...`,
  ## or transforming model parameters
  ##############################
  
  ## Somewhat redundant example
  parameter_names <- parTab$names
  
  ##############################
  ## This is where you would put your own model code
  ##############################
  f <- function(pars){
    names(pars) <- parameter_names
    sim_length <- pars["sim_length"]
    starting_EIR <- pars["starting_EIR"]
    
    
    simparams <- get_parameters(
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
    
    simparams <- set_equilibrium(simparams, starting_EIR)
    
    bednetparams <- simparams
    
    bednet_events = data.frame(
      timestep = 90,
      name=c("Bednets")
    )
    
    bednetparams <- set_bednets(
      bednetparams,
      timesteps = bednet_events$timestep,
      coverages = c(.8),
      retention = 5 * year,
      dn0 = matrix(c(.533), nrow=1, ncol=1),
      rn = matrix(c(.56), nrow=1, ncol=1),
      rnm = matrix(c(.24), nrow=1, ncol=1),
      gamman = 2.64 * 365
    )
    
    print("Sim starting")
    output <- run_simulation(sim_length, bednetparams)
    print("Sim finished")
    
    model_prev <- (output$n_detect_730_3650 / output$n_730_3650)[85]
    
    ## Note the use of closures to capture the `data` argument
    lik <- dpois(round(data * output$n_730_3650[85]), output$n_detect_730_3650[85], log = TRUE)
    if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
    lik
  }
  return(f)
}

## Note that we've given it a prior function too. lazymcmc will 
## combine this function pointer with the likelihood function.
my_prior <- function(pars){
  a <- dnorm(pars["starting_EIR"],50,10, log = TRUE)
  return(a)
}

# Plotting functions
plot_prevalence <- function(output) {
  ggplot(output) + geom_line(
    aes(x = timestep, y = (n_detect_730_3650 / n_730_3650))) + 
    labs(x = "timestep", y = "PfPR2-10")
}

add_intervention_lines <- function(plot, events) {
  plot + geom_vline(
    data = events,
    mapping = aes(xintercept=timestep),
    color="blue"
  ) + geom_text(
    data = events,
    mapping = aes(x = timestep, y = 0, label = name),
    size = 4,
    angle = 90,
    vjust = -0.4,
    hjust = 0
  )
}