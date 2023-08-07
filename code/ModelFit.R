# This script ...

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import data, create df with abundances----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

#---- 1.1. Import packages ----
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(rstan)
#install.packages("loo")
library(loo) # Efficient approximate leave-one-out cross-validation
#install.packages("HDInterval")
library(HDInterval)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("dplyr")
library(dplyr)
library(ggpubr)
library(ggplot2)
library(odeintr)
library(brms)
library(bbmle)
library(reshape2)
#rstan_options(auto_write = TRUE)

#---- 1.2. Import the Data ----
#####################################
Ricker = F
Stouffer= T
Malyon= F

# define initial conditions and model parameters
# define all parameter values
# paired values are always ordered (i, j)

run = T
for( scenario in c("low","medium","high")){

# initial conditions for viable seeds in seed bank
N0 <- c(200,10)
time.exp.max <- 8 # number of years of data collected
time.exp.min <- 5
# number of years to simulate
nyears <- 100
params.low <- list(
  alpha_ij = 1 , #competitive effect of j on i
  alpha_ji = 0.50 #competitive effect of i on j
)

params.medium <- list(
  alpha_ij = 0.50 , #competitive effect of j on i
  alpha_ji = 0.50 #competitive effect of i on j

)
params.high <- list(
  alpha_ij = 0.5 , #competitive effect of j on i
  alpha_ji = 1 #competitive effect of i on j
)
params <- append(get(paste0("params.",scenario)),
                 list( # constant parameters
                   mu    = c(0.10, 0.10), # mortality rate of seeds
                   nu    = c(0.10, 0.10), # mortality rate of ind
                   T     = 0.50,
                   gamma = c(0.1, 0.1), # germination rate of seeds
                   r     = c(2.00, 2.00), # intrinsic growth rate
                   K     = c(25.0, 25.0), # carrying capacity
                   beta  = c(0.2, 0.2), # biomass of germinant 
                   phi   = c(5,5) # conversion rate from biomass to seed
                 ))

params.plant <- params #for plant growth phase

params.seed <- append(get(paste0("params.",scenario)), # for seed growth phase
                      list( # constant parameters
                        mu    = c(0.10, 0.10), # mortality rate of seeds
                        nu    = c(0.10, 0.10), # mortality rate of ind
                        T     = 0.50,
                        gamma = c(0.7, 0.7), # germination rate of seeds
                        r     = c(0.00, 0.00), # intrinsic growth rate
                        K     = c(25.0, 25.0), # carrying capacity
                        beta  = c(0.2, 0.2) # biomass of germinant 
                      ))

#for seed germination phase
source("code/GenerateSimData_wrapper.R")

#Generate.experimental.outcomes <- read.csv("results/Generate.experimental.outcomes.csv")
#Generate.simulated.data <- read.csv("results/Generate.simulated.data.csv")

#---- 1.3. Set parameters ----
simdata <- simulated.data
Alphadistribution.neighbours <- data.frame()
Fecunditydistribution <- data.frame()
PostFecunditydistribution <- data.frame()
nspec <- length(c("i","j"))
for(Code.focal in c("i","j")){ #,"j"
  for (function.int in c(1:4)){ # c(1:4)
    
    print(paste(scenario, Code.focal,", function",function.int))
    
    function.vec <- c(0,0,0,0)
    function.vec[function.int] <- 1

# subset for one Code.focal species 
    
SpDataFocal <- simdata[which(simdata$focal == Code.focal),]
#SpDataFocal <- SpDataFocal[which(simulated.data$time>=time.exp.min & 
#                                   simulated.data$time < time.exp.max),]

#SpDataFocal <- simdata
#SpDataFocal <- SpDataFocal[which(SpDataFocal$focal == Code.focal),]

SpDataFocal <- SpDataFocal[complete.cases(SpDataFocal$fecundity),] 
#SpDataFocal$seeds[is.na(SpDataFocal$seeds)] <- 0

# Next continue to extract the data needed to run the model. 
N <- as.integer(nrow(SpDataFocal))
Fecundity <- as.integer(SpDataFocal$fecundity)  

#---- 2. ABUDANCE MATRIX----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2.1. Interaction (direct) matrix of plant with COMP ----

# Now calculate the total number of plant species to use for the model, discounting
#       any species columns with 0 abundance. Save a vector of the species names
#       corresponding to each column for easy matching later.
AllSpNames <- names(SpDataFocal)[!names(SpDataFocal) %in% c("time","focal","seeds.i",
                                                            "seeds.j","fecundity",
                                                            "seeds","biomass")]

AllSpAbunds <- SpDataFocal %>% 
  dplyr::select(all_of(AllSpNames))

SpTotals <- colSums(AllSpAbunds)
SpToKeep <- SpTotals > 0

S <- sum(SpToKeep)
SpMatrix <- matrix(NA, nrow = N, ncol = S)
i <- 1
for(s in 1:ncol(AllSpAbunds)){
  if(SpToKeep[s] == 1){
    SpMatrix[,i] <- AllSpAbunds[,s]
    i <- i + 1
  }else{next}
}
#SpMatrix <-round((SpMatrix/max(SpMatrix))*100) #scale all the interaction between 0 and 100
#if(max(SpMatrix) == 100){print("scale SpMatrix_plant correct")}

SpNames <- AllSpNames[SpToKeep]

#assign(paste0("SpNames_",FocalPrefix),
#     SpNames)
Intra <- ifelse(SpNames == paste0("plants.",Code.focal), 1, 0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 3. BAYES FIT----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- 3.1. Set up summary interactions df and parameters ---- 

run_estimation <- 1
alphaFunct1 <- function.vec[1]
alphaFunct2 <- function.vec[2]
alphaFunct3 <- function.vec[3]
alphaFunct4 <- function.vec[4]

alpha.function <- paste0("function_",which(function.vec==1))

Nmax <- c(SpDataFocal[which.max(SpDataFocal$fecundity),"plants.i"], 
          SpDataFocal[which.max(SpDataFocal$fecundity),"plants.j"])

DataVec <- c("N", "S","Nmax",
             "Fecundity", "SpMatrix",
             "Intra","run_estimation","alphaFunct1",
             "alphaFunct2","alphaFunct3","alphaFunct4")

##---- 3.2. Run  final fit ----
# Now run a fianl fit of the model to assess parameter 
print("Final Fit beginning")

#install.packages("codetools")
library("codetools")
options(mc.cores = parallel::detectCores())
list.init <- function(...)list(lambdas= array(abs(as.numeric(rnorm(1,mean = 2,
                                                                   sd= 1), dim = 1))))
if(run == T){
FinalFit <- rstan::stan(file = "code/DensityFunct_BH_Final.stan", 
                  data = DataVec,
                  init=list.init,
                  control=list(max_treedepth=10),
                  warmup= 1500,
                  iter = 3000, 
                  chains = 4,
                  seed=1616)

# for function 3 and 4 with divergence, rerun by giving my informative init.value from function 2 values
if( alpha.function >=3 & max(summary(FinalFit)$summary[,"Rhat"],na.rm =T) >1.1 ){
  load(paste0("results/FinalFit_",scenario,"_",Code.focal,"_function_",2,".rds"))
  lambdas = array(c(summary(FinalFit)$summary['lambdas[1]',"mean"]), dim = 1)
  alpha_initial = array(c(summary(FinalFit)$summary['alpha_initial[1]',"mean"],
                          summary(FinalFit)$summary['alpha_initial[2]',"mean"]), dim = nspec)
  alpha_slope = array(c(summary(FinalFit)$summary['alpha_slope[1]',"mean"],
                        summary(FinalFit)$summary['alpha_slope[2]',"mean"]), dim = nspec)
  
  list.init <- function(...)list(lambdas= lambdas,
                                 alpha_initial = alpha_initial,
                                 alpha_slope =  alpha_slope )
  remove(FinalFit)
  FinalFit <- rstan::stan(file = "code/DensityFunct_BH_Final.stan", 
                          data = DataVec,
                          init=list.init,
                          control=list(max_treedepth=12),
                          warmup= 1500,
                          iter = 3000, 
                          chains = 4,
                          seed=16)
  remove(lambdas,alpha_initial,alpha_slope)
}

save(file= paste0("results/FinalFit_",scenario,"_",Code.focal,"_",alpha.function,".rds"),
     FinalFit)
}
load(paste0("results/FinalFit_",scenario,"_",Code.focal,"_",alpha.function,".rds"))

FinalPosteriors <- rstan::extract(FinalFit)

print("Final Fit done")

#---- 3.3. Final fit posterior check and behavior checks---- 

##### Diagnostic plots and post prediction 
pdf(paste0("figures/FinalFit_",scenario,"_",Code.focal,"_",alpha.function,".pdf"))
# Extract pointwise log-likelihood
# using merge_chains=FALSE returns an array, which is easier to 
# use with relative_eff()

# Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
source("code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
# check the distribution of Rhats and effective sample sizes 
##### Posterior check
stan_post_pred_check(FinalPosteriors,"F_hat",Fecundity,
                     paste0("results/PostFec_",scenario,"_",
                            Code.focal,"_",alpha.function,".csv.gz")) 

log_lik_2 <- loo::extract_log_lik(FinalFit, 
                                  parameter_name = "F_sim", 
                                  merge_chains = F)

r_eff <- loo::relative_eff(exp(log_lik_2), cores = 2) 

loo_1 <- loo::loo(log_lik_2 , r_eff = r_eff, cores = 2)

print(loo_1)

# N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
hist(summary(FinalFit)$summary[,"Rhat"],
     main = paste("Finat Fit: Histogram of Rhat for",
                  Code.focal," and ",alpha.function))
hist(summary(FinalFit)$summary[,"n_eff"],
     main = paste("Finat Fit: Histogram of Neff for",
                  Code.focal," and ",alpha.function))

# plot the corresponding graphs
#stan_model_check(FinalFit,
 #                param =c('lambdas','c','alpha_initial','alpha_slope','c'))
# Next check the correlation among key model parameters and identify any
#pairs(FinalFit, pars = c("lambdas",'alpha_initial','alpha_slope','c'))

# functions from Rstan pacakges
trace <- stan_trace(FinalFit, pars=c('lambdas','c','alpha_initial','alpha_slope','c'),
           inc_warmup = TRUE)
print(trace)
dens <- stan_dens(FinalFit, 
          pars=c('lambdas','c','alpha_initial','alpha_slope','c'))
print(dens)
splot <- stan_plot(FinalFit, 
          pars=c('lambdas','c','alpha_initial','alpha_slope','c'))
print(splot)
sampler_params <- get_sampler_params(FinalFit, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)
pairs(FinalFit, pars = c("lambdas",'alpha_initial',
                         'alpha_slope','c'))
#The below-diagonal intersection and the above-diagonal intersection of the same two variables should have distributions that are mirror images of each other.
#Any yellow points would indicate transitions where the maximum treedepth__ was hit, and red points indicate a divergent transition.


#print(loo_1)

dev.off()

library(car)

#---- 3.4. Extraction interactions coefficients---

density.comp <- data.frame(observations= c(1:nrow(SpDataFocal)),
                           species.i = SpDataFocal$plants.i,
                           species.j = SpDataFocal$plants.j)

Alphadistribution.i <- tibble()
Alphadistribution.j <- tibble()

Alphadistribution.i <- data.frame(FinalPosteriors$alpha_value[,,1])
names(Alphadistribution.i) <- c(1:nrow(SpDataFocal))
Alphadistribution.i <- gather(Alphadistribution.i, key="observations",value="alpha.i")
Alphadistribution.i$observations <- as.numeric(Alphadistribution.i$observations)
Alphadistribution.i <- full_join(Alphadistribution.i,density.comp, 
                                 by=c("observations"))

Alphadistribution.j <- as.data.frame(FinalPosteriors$alpha_value[,,2])
names(Alphadistribution.j) <- c(1:nrow(SpDataFocal))
Alphadistribution.j <- gather(Alphadistribution.j, key="observations",value="alpha.j")
Alphadistribution.j$observations <- as.numeric(Alphadistribution.j$observations)
Alphadistribution.j <- full_join(Alphadistribution.j,density.comp, by=c("observations"))


Alphadistribution.i <- Alphadistribution.i %>% 
  group_by(species.i) %>% summarise_at("alpha.i",  list(mean = mean, sd = sd))
                                                           
Alphadistribution.i <- data.frame(abundance.neighbours = Alphadistribution.i$species.i,
                                         alpha_mean = Alphadistribution.i$mean,
                                         alpha_sd= Alphadistribution.i$sd,
                                         neighbours= "species i",focal = paste("species",Code.focal),
                                         density.function = alpha.function,
                                         Nmax = Nmax[1])


Alphadistribution.j <- Alphadistribution.j %>% 
  group_by(species.j) %>% summarise_at("alpha.j",  list(mean = mean, sd = sd))

Alphadistribution.j <- data.frame(abundance.neighbours = Alphadistribution.j$species.j,
                                  alpha_mean = Alphadistribution.j$mean,
                                  alpha_sd= Alphadistribution.j$sd,
                                  neighbours= "species j",focal = paste("species",Code.focal),
                                  density.function = alpha.function,
                                  Nmax = Nmax[2])


Alphadistribution.neighbours <- bind_rows(Alphadistribution.neighbours,Alphadistribution.i,Alphadistribution.j)

#---- 3.4. Extraction fecundity---

Fecunditydistribution.n <- FinalPosteriors %>% 
  as.data.frame() %>% 
  dplyr::select(contains("F_sim")) %>%
  gather( key="obervation",value="seed")
Fecunditydistribution.n$focal = paste("species",Code.focal)
Fecunditydistribution.n$density.function = alpha.function

Fecunditydistribution <- bind_rows(Fecunditydistribution,Fecunditydistribution.n)

PostFecunditydistribution.n  <- read.csv(paste0("results/PostFec_",scenario,"_",Code.focal,"_",alpha.function,".csv.gz"))
PostFecunditydistribution.n$focal <- Code.focal
PostFecunditydistribution.n$alpha.function <- alpha.function

PostFecunditydistribution <- bind_rows(PostFecunditydistribution,PostFecunditydistribution.n)



  }
}

save(Alphadistribution.neighbours, file = paste0("results/Alphadistribution.",scenario,".neighbours.csv.gz"))
save(Fecunditydistribution, file = paste0("results/Fecunditydistribution.",scenario,".csv.gz"))
save(PostFecunditydistribution , file = paste0("results/PostFecunditydistribution.",scenario,".csv.gz"))

}





