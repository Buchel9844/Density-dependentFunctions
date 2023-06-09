
# libraries required to run the code
library(odeintr)

# to use/run the 2sp model
#source('lib/model.2sp.R')

#################################################
# define initial conditions and model parameters
#################################################

# initial conditions for viable seeds in seed bank
N0 <- c(770, 25000)

# number of years to simulate
nyears <- 10

# define all parameter values for seed germination phase
# paired values are always ordered (i, j)
params.seed <- list(
  T     = 0.50,
  gamma = c(0.10, 0.01),# germination rate of seeds
  mu    = c(0.10, 0.20), # mortality rate of seeds
  nu    = c(0.00, 0.00), # mortality rate of ind
  r     = c(0.00, 0.00), # intrinsic growth rate
  K     = c(100.0, 250.0), # carrying capacity
  beta  = c(0.02, 0.02), # biomass of germinant 
  alpha_ij = -0.1,
  alpha_ji = 0.50
) 

# define all parameter values for plant growth phase
# paired values are always ordered (i, j)
params.plant <- list(
  T     = 0.50,
    gamma = c(0.00, 0.001), # germination rate of seeds
    mu    = c(0.10, 0.20), # mortality rate of seeds
    nu    = c(0.20, 0.05), # mortality rate of ind
    r     = c(10.00, 10.00), # intrinsic growth rate
    K     = c(100.0, 250.0), # carrying capacity
    beta  = c(0.2, 0.2), # biomass of germinant 
    alpha_ij = 0.05, #competitive effect of j on i
    alpha_ji = 0.50, #competitive effect of i on j
    phi   = c(10,25) # conversion rate from biomass to seed
  )

################################################
# simulate continuous-time population dynamics
# using above initial conditions and parameters
################################################

# initial conditions including viable seeds, plants, and biomass
population.dynamics <- matrix(
  c(
    0,
    N0[1], 0, 0,
    N0[2], 0, 0
  ),
  byrow=TRUE,
  nrow=1
)
colnames(population.dynamics) <- c(
  "Time",
  "Seeds i",
  "Plants i",
  "Biomass i",
  "Seeds j",
  "Plants j",
  "Biomass j"
)

# run simulation for the specified number of years
for(year in 0:(nyears-1)){
  # set the parameters for phase 1
  fecundity_dynamics_set_params(
    gamma_i = params.seed$gamma[1],
    mu_i = params.seed$mu[1],
    nu_i = params.seed$nu[1],
    r_i = params.seed$r[1],
    K_i = params.seed$K[1],
    beta_i = params.seed$beta[1],
    gamma_j = params.seed$gamma[2],
    mu_j = params.seed$mu[2],
    nu_j = params.seed$nu[2],
    r_j = params.seed$r[2],
    K_j = params.seed$K[2],
    beta_j = params.seed$beta[2],
    alpha_ij = params.seed$alpha_ij,
    alpha_ji = params.seed$alpha_ji
  )
  # simulate phase 1
  NPBseed <- fecundity_dynamics(
    init=population.dynamics[nrow(population.dynamics),2:7],
    duration=params.seed$T[1],
    step_size=params.seed$T[1]/1000.,
    start=year
  )
  # add phase 1 dynamics to pop dyn container
  population.dynamics <- rbind(population.dynamics, as.matrix(NPBseed))
  
  # set the parameters for phase 2
  fecundity_dynamics_set_params(
    gamma_i = params.plant$gamma[1],
    mu_i = params.plant$mu[1],
    nu_i = params.plant$nu[1],
    r_i = params.plant$r[1],
    K_i = params.plant$K[1],
    beta_i = params.plant$beta[1],
    gamma_j = params.plant$gamma[2],
    mu_j = params.plant$mu[2],
    nu_j = params.plant$nu[2],
    r_j = params.plant$r[2],
    K_j = params.plant$K[2],
    beta_j = params.plant$beta[2],
    alpha_ij = params.plant$alpha_ij,
    alpha_ji = params.plant$alpha_ji
  )
  # simulate phase 2
  NPBplant <- fecundity_dynamics(
    init=population.dynamics[nrow(population.dynamics),2:7],
    duration=params.plant$T[1],
    step_size=params.plant$T[1]/1000.,
    start=population.dynamics[nrow(population.dynamics),1]
  )
  # add phase 2 dynamics to pop dyn container
  population.dynamics <- rbind(population.dynamics, as.matrix(NPBplant))
  
  # convert biomass to seeds, kill all plants, remove their biomass
  end.state <- population.dynamics[nrow(population.dynamics),]
  end.state["Seeds i"] <- end.state["Seeds i"] + params.plant$phi[1] * end.state["Biomass i"]
  end.state["Seeds j"] <- end.state["Seeds j"] + params.plant$phi[2] * end.state["Biomass j"]
  end.state["Plants i"] <- end.state["Plants j"] <- 0
  end.state["Biomass i"] <- end.state["Biomass j"] <- 0
  
  # add this last step to the population dynamics
  population.dynamics <- rbind(population.dynamics, end.state)
}

###################################################
# plot two-species population-dynamics time series 
###################################################

dev.new(width=16,height=7)
par(oma=c(0,0,0,0), mar=c(4.5,5,4,2))
layout(mat = matrix(
  1:6,
  byrow=TRUE,
  nrow = 2,
  ncol = 3
),
heights = c(6,6),
widths = rep(3,6)
)

par(cex.lab=1.5,cex.axis=1.2)
par(lwd=1.5)

# seeds of i
plot(
  population.dynamics[,"Time"],
  population.dynamics[,"Seeds i"],
  type='l',
  col=1,
  lty=1,
  xlab="Time",
  ylab=expression("Viable seeds in the seed bank ("*italic(N[i])*")"),
  xlim=c(0,nyears)
)

# plants of i
plot(
  population.dynamics[,"Time"],
  population.dynamics[,"Plants i"],
  type='l',
  col=1,
  lty=1,
  xlab="Time",
  ylab=expression("Plants ("*italic(p[i])*")"),
  xlim=c(0,nyears)
)

# biomass of i
plot(
  population.dynamics[,"Time"],
  population.dynamics[,"Biomass i"],
  type='l',
  col=1,
  lty=1,
  xlab="Time",
  ylab=expression("Plant biomass ("*italic(B[i])*")"),
  xlim=c(0,nyears)
)

# seeds of j
plot(
  population.dynamics[,"Time"],
  population.dynamics[,"Seeds j"],
  type='l',
  col=2,
  lty=2,
  xlab="Time",
  ylab=expression("Viable seeds in the seed bank ("*italic(N[j])*")"),
  xlim=c(0,nyears)
)

# plants of j
plot(
  population.dynamics[,"Time"],
  population.dynamics[,"Plants j"],
  type='l',
  col=2,
  lty=2,
  xlab="Time",
  ylab=expression("Plants ("*italic(p[j])*")"),
  xlim=c(0,nyears)
)

# biomass of j
plot(
  population.dynamics[,"Time"],
  population.dynamics[,"Biomass j"],
  type='l',
  col=2,
  lty=2,
  xlab="Time",
  ylab=expression("Plant biomass ("*italic(B[j])*")"),
  xlim=c(0,nyears)
)
