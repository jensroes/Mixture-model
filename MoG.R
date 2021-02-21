# Load packages
library(tidyverse)
library(rstan)
library(MASS) # using mvrnorm
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
n_chain = n_core = 3 # number of cores/chains
iterations = 6e3
set.seed(1)

# Get some fake data (see script for parameter values)
source("MoG_data.R"); data

# Data as list
dat <- 
  within(list(), {
    N <- nrow(data)
    y <- data$value
    condition <- data$condition # or whatever grouping variable
    subj <- as.integer(data$id)
    S <- max(as.integer(data$id))
#    items <- as.integer(data$item)
#    I <- max( as.integer(data$item))
  }  ); str(dat)

# Initialise start values
start <-  function(chain_id = 1){
    list(   beta = 6
          , delta =.1
          , theta = rep(0,2)
          , sigma = 1
          , sigma_diff = .1
          , u = rep(0, dat$S)
          , sigma_u = 0.1 ) }

start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )

# Check compiling
mog <- stan(file = "MoG.stan", data=dat, chains=0)

# Fit model
m <- stan(fit = mog, data = dat, init = start_ll,
         iter = iterations, warmup= iterations/2,
         chains = n_chain, cores = n_core, 
         refresh = 2000, seed = 365,
         control = list(max_treedepth = 16,
                        adapt_delta = 0.99,
                        stepsize = 2))

# Save posterior samples
#saveRDS(m,
#        file="mog.rda",
#        compress="xz")

#m <- readRDS("mog.rda")

param <- c("beta", "delta", "beta", "prob", "sigma")
traceplot(m, param)
summary(m, param, prob = c(.025, .975))$summary %>% round(2)
