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

Ricker = T
Stouffer= F
Malyon= F
source("code/GenerateSimData_wrapper.R")


#---- 1.3. Set parameters ----
names(simdata) <- c("focal","fecundity","i","j")
Alphadistribution.neighbours.Ricker <- data.frame()
Fecunditydistribution.Ricker <- data.frame()
PostFecunditydistribution.Ricker <- data.frame()
for(Code.focal in c("i","j")){ #,"j"
  for (function.int in c(1:4)){ # c(1:4)
    print(paste(Code.focal,", function",function.int))
    
    function.vec <- c(0,0,0,0)
    function.vec[function.int] <- 1
    
    # subset for one Code.focal species 
    SpDataFocal <- simdata[which(simdata$focal == Code.focal),]
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
    FinalFit <- rstan::stan(file = "code/DensityFunct_BH_Final.stan", 
                            data = DataVec,
                            init="random",
                            warmup= 500,
                            iter = 1000, 
                            chains = 3)
    
    
    save(file= paste0("results/FinalFit_Ricker_",Code.focal,"_",alpha.function,".rds"),
         FinalFit)
    
    load(paste0("results/FinalFit_Ricker_",Code.focal,"_",alpha.function,".rds"))
    
    FinalPosteriors <- rstan::extract(FinalFit)
    
    print("Final Fit done")
    
    #---- 3.3. Final fit posterior check and behavior checks---- 
    
    ##### Diagnostic plots and post prediction 
    pdf(paste0("figures/FinalFit_Ricker_",Code.focal,"_",alpha.function,".pdf"))
    # Extract pointwise log-likelihood
    # using merge_chains=FALSE returns an array, which is easier to 
    # use with relative_eff()
    
    # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
    source("code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
    # check the distribution of Rhats and effective sample sizes 
    ##### Posterior check
    stan_post_pred_check(FinalPosteriors,"F_hat",Fecundity,
                         paste0("results/PostFec_Ricker_",Code.focal,"_",alpha.function,".csv.gz")) 
    
    # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
    hist(summary(FinalFit)$summary[,"Rhat"],
         main = paste("Finat Fit: Histogram of Rhat for",
                      Code.focal," and ",alpha.function),
         xlab = "Rhat")
    hist(summary(FinalFit)$summary[,"n_eff"],
         main = paste("Finat Fit: Histogram of Neff for",
                      Code.focal," and ",alpha.function),
         xlab = "Rhat")
    
    # plot the corresponding graphs
    #stan_model_check(FinalFit,
    #                param =c('lambdas','c','alpha_initial','alpha_slope','c'))
    # Next check the correlation among key model parameters and identify any
    #pairs(FinalFit, pars = c("lambdas",'alpha_initial','alpha_slope','c'))
    
    # functions from Rstan pacakges
    stan_trace(FinalFit, pars=c('lambdas','c','alpha_initial','alpha_slope','c'),
               inc_warmup = TRUE)
    stan_dens(FinalFit, pars=c('lambdas','c','alpha_initial','alpha_slope','c'))
    stan_plot(FinalFit, pars=c('lambdas','c','alpha_initial','alpha_slope','c'))
    
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
                               species.i = SpDataFocal$i,
                               species.j = SpDataFocal$j)
    
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
                                      density.function = alpha.function)
    
    
    Alphadistribution.j <- Alphadistribution.j %>% 
      group_by(species.j) %>% summarise_at("alpha.j",  list(mean = mean, sd = sd))
    
    Alphadistribution.j <- data.frame(abundance.neighbours = Alphadistribution.j$species.j,
                                      alpha_mean = Alphadistribution.j$mean,
                                      alpha_sd= Alphadistribution.j$sd,
                                      neighbours= "species j",focal = paste("species",Code.focal),
                                      density.function = alpha.function)
    
    
    Alphadistribution.neighbours.Ricker <- bind_rows(Alphadistribution.neighbours.Ricker,Alphadistribution.i,Alphadistribution.j)
    
    #---- 3.4. Extraction fecundity---
    
    Fecunditydistribution.Ricker.n <- FinalPosteriors %>% 
      as.data.frame() %>% 
      dplyr::select(contains("F_sim")) %>%
      gather( key="obervation",value="seed")
    Fecunditydistribution.Ricker.n$focal = paste("species",Code.focal)
    Fecunditydistribution.Ricker.n$density.function = alpha.function
    
    Fecunditydistribution.Ricker <- bind_rows(Fecunditydistribution.Ricker,Fecunditydistribution.Ricker.n)
    
    PostFecunditydistribution.Ricker.n  <- read.csv(paste0("results/PostFec_Ricker_",Code.focal,"_",alpha.function,".csv.gz"))
    PostFecunditydistribution.Ricker.n$focal <- Code.focal
    PostFecunditydistribution.Ricker.n$alpha.function <- alpha.function
    
    PostFecunditydistribution.Ricker <- bind_rows(PostFecunditydistribution.Ricker,PostFecunditydistribution.Ricker.n)
    
    
    
  }
}

save(Alphadistribution.neighbours.Ricker, file = "results/Alphadistribution.neighbours.Ricker.csv.gz")
save(Fecunditydistribution.Ricker, file = "results/Fecunditydistribution.Ricker.csv.gz")
save(PostFecunditydistribution.Ricker , file = "results/PostFecunditydistribution.Ricker.csv.gz")







