# This script will run the empirical model fits for each focal species and each
#       environmental covariate. A separate script will then make the empirical
#       figures for the manuscript

# https://mc-stan.org/docs/2_29/functions-reference/index.html#overview
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import data, create df with abundances----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

#---- 1.1. Import packages ----
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(rstan)
#install.packages("HDInterval")
library("HDInterval")
#install.packages("tidyverse")
library("tidyverse")
#install.packages("dplyr")
library(dplyr)
library(ggpubr)
library(ggplot2)
#rstan_options(auto_write = TRUE)

#---- 1.2. Import the Data ----
setwd("/Users/lisabuche/Projects/Density-dependentFunctions")
#### Stouffer's model
set.seed(1236) #to create reproducible results
source("code/GenerateSimData-Stouffer.R") # Generate Simulated data in "simulated.data" df and 
write_csv(simulated.data, 
           file = "results/simulated.data.csv")
write_csv(experimental.outcomes, 
           file = "results/experimental.outcomes.csv")
head(simulated.data)
ggsave("figures/simulated.seed.density.pdf",
plot = ggplot(simulated.data, aes(x=fecundity, fill=focal, color=focal)) + geom_density(alpha=0.5) + 
  scale_x_continuous(trans='log2',
                     name = "Viable seed, fecundity per individuals (ln transformed)") + 
  theme_bw() 
)

biomass.sim <- c("focal species i","focal species j")
names(biomass.sim) <- c("biomass.i", "biomass.j")
ggsave("figures/simulated.biomass.abundance.pdf",
       plot = ggplot(gather(gather(experimental.outcomes, biomass.i, biomass.j,
                            key="focal.biomass", value="biomass"),
                            ending.plants.i, ending.plants.j,
                            key="focal.ending.plants", value="ending.plants")) +
         geom_smooth(aes(y=biomass, x=per.capita.fecundity.i,
                         color=focal.ending.plants,
                         fill=focal.ending.plants)) +
         facet_grid(.~focal.biomass, 
                    labeller = labeller(focal.biomass = biomass.sim))+ 
         ylim(c(0,1)) +theme_bw() + guides(fill="none") +
         scale_color_discrete(name = "Neighbour \nbiomass",labels=c("species i", "species j")) +
         xlab("abundance focal species") +
         labs(caption = expression("with interaction coefficients"~alpha[ij]==0.9~"and"~alpha[ji]==0.7)) +
         theme(plot.caption.position = "plot",
               plot.caption = element_text(hjust = 0,size=11))
    
)

ggplot(experimental.outcomes) +
  geom_smooth(aes(y=per.capita.fecundity.j, x=per.capita.fecundity.i))
ggplot(experimental.outcomes) +
  geom_point(aes(y=biomass.j, x=per.capita.fecundity.i)) 
#### Ricker model Data
source("code/simul_data.R")
set.seed(16)
simul_data_ricker <- simul_data(S = 2,   # number of focal groups / species
  K = 2,   # number of neighbour focals 
  pI = 0.1) # Generate Simulated data in "simulated.data" df and
simdata <- simul_data_ricker$simdata
cols.num <- c("seeds","i","j")
simdata[cols.num] <- sapply(simdata[cols.num],as.numeric)

ggsave("figures/Ricker.fecundity.pdf",
       plot = ggplot(simdata, aes(x=log(seeds), color=as.factor(focal))) + 
         geom_density(alpha=0.5) + 
         scale_x_continuous(name = "Viable seed, fecundity per individuals (log transformed)") + 
         theme_bw() 
)
simul_data_ricker$sim_alpha_specific
#---- 1.3. Set parameters ----
Alphadistribution.neighbours <- data.frame()
Fecunditydistribution <- data.frame()
PostFecunditydistribution <- data.frame()
for(Code.focal in c("i","j")){
  for (function.int in c(2:4)){
    print(paste(Code.focal,", function",function.int))
    
    function.vec <- c(0,0,0,0)
    function.vec[function.int] <- 1

# subset for one Code.focal species 
SpDataFocal <- simulated.data[which(simulated.data$focal == Code.focal),]
#SpDataFocal <- simdata
SpDataFocal <- SpDataFocal[which(SpDataFocal$focal == Code.focal),]

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
                                                            "seeds.j","fecundity","seeds","biomass")]

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
Intra <- ifelse(SpNames == Code.focal, 1, 0)

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

DataVec <- c("N", "S",
             "Fecundity", "SpMatrix",
             "Intra","run_estimation","alphaFunct1",
             "alphaFunct2","alphaFunct3","alphaFunct4")

##---- 3.2. Run  final fit ----
# Now run a fianl fit of the model to assess parameter 
print("Final Fit beginning")

#install.packages("codetools")
library("codetools")
options(mc.cores = parallel::detectCores())
FinalFit <- stan(file = "code/DensityFunct_BH_Final.stan", 
                 data = DataVec,
                  init="random",
                  warmup= 500,
                  iter = 1000, 
                  chains = 3)


save(file= paste0("results/FinalFit_",Code.focal,"_",alpha.function,".rds"),
     FinalFit)

load(paste0("results/FinalFit_",Code.focal,"_",alpha.function,".rds"))

FinalPosteriors <- rstan::extract(FinalFit)

print("Final Fit done")

#---- 3.3. Final fit posterior check and behavior checks---- 

##### Diagnostic plots and post prediction 
pdf(paste0("figures/FinalFit_",Code.focal,"_",alpha.function,".pdf"))
# Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
source("code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
# check the distribution of Rhats and effective sample sizes 
##### Posterior check
stan_post_pred_check(FinalPosteriors,"F_hat",Fecundity,
                     paste0("results/PostFec_",Code.focal,"_",alpha.function,".csv.gz")) 


# N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
hist(summary(FinalFit)$summary[,"Rhat"],
     main = paste("Finat Fit: Histogram of Rhat for",
                  Code.focal," and ",alpha.function))
hist(summary(FinalFit)$summary[,"n_eff"],
     main = paste("Finat Fit: Histogram of Neff for",
                  Code.focal," and ",alpha.function))

# plot the corresponding graphs
stan_model_check(FinalFit,
                 param =c('lambdas','c','alpha_initial','alpha_slope','c'))

# Next check the correlation among key model parameters and identify any
pairs(FinalFit, pars = c("lambdas",'alpha_initial','alpha_slope','c'))


dev.off()

#---- 3.4. Extraction interactions coefficients---

density.comp <- data.frame(observations= c(1:nrow(SpDataFocal)),
                           species.i = SpDataFocal$plants.i,
                           species.j = SpDataFocal$plants.j)

Alphadistribution.i <- tibble()
Alphadistribution.j <- tibble()

Alphadistribution.i <- data.frame(FinalPosteriors$alpha_function_eij[,,1])
names(Alphadistribution.i) <- c(1:nrow(SpDataFocal))
Alphadistribution.i <- gather(Alphadistribution.i, key="observations",value="alpha.i")
Alphadistribution.i$observations <- as.numeric(Alphadistribution.i$observations)
Alphadistribution.i <- full_join(Alphadistribution.i,density.comp, 
                                 by=c("observations"))

Alphadistribution.j <- as.data.frame(FinalPosteriors$alpha_function_eij[,,2])
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
                                         density.function = alpha.function)


Alphadistribution.j <- Alphadistribution.j %>% 
  group_by(species.j) %>% summarise_at("alpha.j",  list(mean = mean, sd = sd))

Alphadistribution.j <- data.frame(abundance.neighbours = Alphadistribution.j$species.j,
                                  alpha_mean = Alphadistribution.j$mean,
                                  alpha_sd= Alphadistribution.j$sd,
                                  neighbours= "species j",focal = paste("species",Code.focal),
                                  density.function = alpha.function)


Alphadistribution.neighbours <- bind_rows(Alphadistribution.neighbours,Alphadistribution.i,Alphadistribution.j)

#---- 3.4. Extraction fecundity---

Fecunditydistribution.n <- FinalPosteriors %>% 
  as.data.frame() %>% 
  dplyr::select(contains("F_sim")) %>%
  gather( key="obervation",value="seed")
Fecunditydistribution.n$focal = paste("species",Code.focal)
Fecunditydistribution.n$density.function = alpha.function

Fecunditydistribution <- bind_rows(Fecunditydistribution,Fecunditydistribution.n)

PostFecunditydistribution.n  <- read.csv(paste0("results/PostFec_",Code.focal,"_",alpha.function,".csv.gz"))
PostFecunditydistribution.n$focal <- Code.focal
PostFecunditydistribution.n$alpha.function <- alpha.function
  
PostFecunditydistribution <- bind_rows(PostFecunditydistribution,PostFecunditydistribution.n)



  }
}

save(Alphadistribution.neighbours, file = "results/Alphadistribution.neighbours.csv.gz")
save(Fecunditydistribution, file = "results/Fecunditydistribution.csv.gz")
save(PostFecunditydistribution , file = "results/PostFecunditydistribution.csv.gz")







