# Code originaly from Daniel Stouffer Github,for his manuscript  the manuscript 
# "A critical examination of models of annual-plant population dynamics and 
# density-dependent fecundity accepted at Methods In Ecology & Evolution by 
# Daniel B. Stouffer (daniel.stouffer@canterbury.ac.nz).
# https://github.com/stoufferlab/annual-plant-dynamics


##################################
# model simulation functions
##################################

# we use this to compile the dynamic model
#install.packages("odeintr")
#install.packages("brms")
library(odeintr)
library(brms)
library(bbmle)
#install.packages("bbmle")


# define the ode in C++ format to use odeintr
# x[0] = seeds in seed bank species i
# x[1] = plants species i
# x[2] = biomass species i
# x[3] = seeds in seed bank species j
# x[4] = plants species j
# x[5] = biomass species j
fecundity.dynamics.sys = '
// seeds in the seed bank species i
dxdt[0] = -gamma_i * x[0] - mu_i * x[0];

// plants that have germinated species i
dxdt[1] = gamma_i * x[0] - nu_i * x[1];

// plant biomass species i
dxdt[2] = x[2] * r_i * (1 - (x[2] + alpha_ij * x[5]) / K_i)
+ beta_i * gamma_i * x[0]
- nu_i * x[2];

// seeds in the seed bank species j
dxdt[3] = -gamma_j * x[3] - mu_j * x[3];

// plants that have germinated species j
dxdt[4] = gamma_j * x[3] - nu_j * x[4];

// plant biomass species j
dxdt[5] = x[5] * r_j * (1 - (alpha_ji * x[2] + x[5]) / K_j)
+ beta_j * gamma_j * x[3] 
- nu_j * x[5];
'

# compile the model into the local environment
odeintr::compile_sys(
  "fecundity_dynamics",
  fecundity.dynamics.sys,
  pars = c(
    "gamma_i",
    "gamma_j",
    "mu_i",
    "mu_j",
    "nu_i",
    "nu_j",
    "r_i",
    "r_j",
    "K_i",
    "K_j",
    "beta_i",
    "beta_j",
    "alpha_ij",
    "alpha_ji"
  ),
  const = TRUE,
  method = 'bsd'
)

#######################################
# response surface experimental design
#######################################

# data frame for model inputs and outputs
experimental.design <- data.frame(
  focal = c(rep("i",5), rep("j", 5)),
  starting.plants.i = round(abs(rnorm(100, mean = 0, sd = 3))), #c(3, 2, 1, 1, 1, 2, 1, 0, 0, 0),#round(abs(rnorm(10, mean = 0, sd = 3))), #BucheL change from a set vector
  starting.plants.j = round(abs(rnorm(100, mean = 0, sd = 3))), #c(3, 2, 1, 1, 1, 2, 1, 0, 0, 0),#round(abs(rnorm(10, mean = 0, sd = 3))), #BucheL change from a set vector
  time = rep(1, 10),
  stringsAsFactors = FALSE
)

# number of experimental replicates to simulate
nreps <- 10

#####################################
# parameters for the fecundity model
#####################################

# define all parameter values
# paired values are always ordered (i, j)
params <- list(
  gamma = c(0.0, 0.0), # germination rate of seeds
  mu    = c(0.0, 0.0), # mortality rate of seeds
  nu    = c(0.0, 0.0), # mortality rate of ind
  r     = c(10.0, 5.0), # intrinsic growth rate
  K     = c(1.0, 1.0), # carrying capacity
  beta  = c(0.001, 0.001),# biomass of germinant 
  phi   = c(500, 250), # convertion rate from biomass to seed
  alpha_ij = 0.50,
  alpha_ji = 0.70
)

########################################################
# simulate the outcomes the response-surface experiment
########################################################

# set the parameters for the model
fecundity_dynamics_set_params(
  gamma_i = params$gamma[1],
  mu_i = params$mu[1],
  nu_i = params$nu[1],
  r_i = params$r[1],
  K_i = params$K[1],
  beta_i = params$beta[1],
  gamma_j = params$gamma[2],
  mu_j = params$mu[2],
  nu_j = params$nu[2],
  r_j = params$r[2],
  K_j = params$K[2],
  beta_j = params$beta[2],
  alpha_ij = params$alpha_ij,
  alpha_ji = params$alpha_ji
)

# add output columns to the experimental design
experimental.outcomes <- experimental.design
newcols <- c(
  "ending.plants",
  "biomass",
  "total.fecundity",
  "per.capita.fecundity"
)
for(cname in newcols){
  for(sp in c("i","j")){
    experimental.outcomes[,paste0(cname,".",sp)] <- NA
  }
}

# run the model for each experimental treatment and save the results
for(rr in 1:nrow(experimental.outcomes)){
  Ni0 <- experimental.outcomes$starting.plants.i[rr]
  Nj0 <- experimental.outcomes$starting.plants.j[rr]
  time <- experimental.outcomes$time[rr]
  
  # starting conditions are seeds in seed bank, plants, and plant biomass
  x0 <- c(
    0,Ni0,Ni0*params$beta[1],
    0,Nj0,Nj0*params$beta[2]
  )
  
  # output order matches the conditions above but after time in first column
  growing.season <- fecundity_dynamics(init=x0, duration=time, step_size=time/1000.)
  
  # name the output columns for convenience
  colnames(growing.season) <- c(
    "Time",
    "seeds.i",
    "plants.i",
    "biomass.i",
    "seeds.j",
    "plants.j",
    "biomass.j"
  )
  
  # use final plants, final biomass, and conversion rate to estimate total and per capita fecundity of i
  plants_i <- growing.season[nrow(growing.season),"plants.i"]
  biomass_i <- growing.season[nrow(growing.season),"biomass.i"]
  seeds_i <- biomass_i * params$phi[1]
  fecundity_i <- seeds_i / plants_i
  
  # use final plants, final biomass, and conversion rate to estimate total and per capita fecundity of j
  plants_j <- growing.season[nrow(growing.season),"plants.j"]
  biomass_j <- growing.season[nrow(growing.season),"biomass.j"]
  seeds_j <- biomass_j * params$phi[2]
  fecundity_j <- seeds_j / plants_j
  
  # fillin the outcomes into the table
  experimental.outcomes[rr,5:12] <- c(
    plants_i,
    plants_j,
    biomass_i,
    biomass_j,
    seeds_i,
    seeds_j,
    fecundity_i,
    fecundity_j
  )
}

##########################################################
# generate random data set based on simulated predictions
##########################################################

# set a random seed for reproducibility
# obtained from random.org (Timestamp: 2022-02-09 20:06:47 UTC)
#set.seed(743761)

# use a poisson model to generate count fecundities based on the above simulations
simulated.data <- do.call(rbind,apply(
  experimental.outcomes,
  1,
  function(x,nreps){
    data.frame(
      focal=rep(x["focal"], nreps),
      plants.i=rep(as.integer(x["starting.plants.i"]), nreps),
      plants.j=rep(as.integer(x["starting.plants.j"]), nreps),
      seeds.i=0,
      seeds.j=0,
      time=as.numeric(x["time"]),
      biomass=as.numeric(x[paste0("biomass.",x["focal"])]),
      fecundity=rpois(nreps, as.numeric(x[paste0("per.capita.fecundity.",x["focal"])])),
      stringsAsFactors=FALSE
    )
  },
  nreps=nreps
))
rownames(simulated.data) <- 1:nrow(simulated.data)

# save the parameters that were used for simulations
simulated.params <- params

