#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import data, create df with competiton and seed distributions----
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
library(tidyr) #fill is part of tidyr
library(lme4)

#---- 1.2. Import the competitive data ----
competition <- read.csv("/Users/lisabuche/Documents/Projects/Facilitation_density/data/Bendering2022_competition.csv",
                        header = T,stringsAsFactors = F, sep=",",
                        na.strings=c("","NA"))

competition <- competition %>%  
  fill(year, day, month, block, unique.plot) %>%   # extend values of year, plot, subplot, etc 
  group_by(year,day,month,block,unique.plot) %>%
  mutate(individual = row_number()) %>%
  ungroup() %>%
  separate_rows(neighbourhood, sep = "/",
                convert = FALSE) %>%  # separate the competitors in ultiple rows
  separate( neighbourhood,sep = "_",
            into = c("competitor", "abundance")) %>%
  spread(key=abundance,value=competitor)

#names(competition) <- tolower(names(competition))

# change the format of the rows to numeric 
cols.num <- names(competition)[which(!names(competition) %in% c("focal"))]
competition <- as.data.frame(competition)
competition[cols.num] <- sapply(competition[cols.num],as.numeric)
# change na values to 0
competition[is.na(competition)] <- 0

#---- 1.2. Import the seeds data ----
fecundity.df <- read.csv("/Users/lisabuche/Documents/Projects/Facilitation_density/data/Bendering2022_seeds.csv",
                         header = T,stringsAsFactors = F, sep=",",
                         na.strings=c("","NA"))

fecundity.df <- fecundity.df %>%  
  fill(year, block, focal)  # extend values of year, plot, subplot, etc 

fecundity.df$block <- as.factor(fecundity.df$block)
seeds.model <- glmer(viable.seeds ~ (1|block) + focal*tot.flowers ,fecundity.df, 
                     family="poisson")
summary(seeds.model) # random block effect not signidicant - not considered further
  
seeds.graph <- ggplot(fecundity.df,aes(y=viable.seeds,x=as.factor(tot.flowers))) + 
    geom_violin() + facet_grid(.~focal,scales="free")
# tot number of flowers on individual not cange per capita fecundity - or not enough data to elaborate that

fecundity.summary <- fecundity.df %>%  
  group_by(year, focal) %>%  
  summarise(mean.seed=mean(viable.seeds), st.dev.seed = sd(viable.seeds),
          n= length(viable.seeds),
          mean.flower=mean(tot.flowers),
          st.dev.seed= sd(tot.flowers))%>%  
  ungroup()


#---- 1.3. join competition with seeds production for each focal  ----

focal.levels <- levels(as.factor(competition$focal))

Spcompetition <- list()

for(focal in focal.levels){
  
  competition.n <- subset(competition, competition$focal == focal) 
  # multiple by a number randomly draw from a normal distribution following seed production of the focal
  competition.n$seeds <-  competition.n$flower*
    rnorm(nrow(competition.n),
          mean=fecundity.summary$mean.seed[which(fecundity.summary$focal==focal)],
          sd = fecundity.summary$st.dev.seed[which(fecundity.summary$focal==focal)])

  #save(competition.n,
  #     paste0("results/Competition_seeds_2022_",focal,".csv"))
  Spcompetition[[focal]] <- competition.n
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2. Run the model for each focal----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(Code.focal in focal.levels){
for (function.int in c(1:4)){
  print(paste(Code.focal,", function",function.int))
  
  function.vec <- c(0,0,0,0)
  function.vec[function.int] <- 1
  
 # data for the focal
  SpDataFocal <- Spcompetition[[focal]]
 
  # Next continue to extract the data needed to run the model. 
  N <- as.integer(nrow(SpDataFocal))
  Fecundity <- as.integer(SpDataFocal$seeds)  
  
  #---- 2. ABUDANCE MATRIX----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #---- 2.1. Interaction (direct) matrix of plant with COMP ----
  
  # Now calculate the total number of plant species to use for the model, discounting
  #       any species columns with 0 abundance. Save a vector of the species names
  #       corresponding to each column for easy matching later.
  AllSpNames <- names(SpDataFocal)[!names(SpDataFocal) %in% c("focal","year","day","month",
                                                              "block","unique.plot",
                                                              "flower","seeds","individual")]
  

  AllSpAbunds <- SpDataFocal %>% 
    dplyr::select(all_of(AllSpNames))%>%
    mutate_at( AllSpNames, as.numeric)
 
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
  print("Natural Fit beginning")
  
  #install.packages("codetools")
  library("codetools")
  options(mc.cores = parallel::detectCores())
  # defiining initial value of lambda to help with computation failure
  list.init <- function(...)list(lambdas = array(as.numeric(rnorm(1,
                                                                  mean=log(fecundity.summary$mean.seed[which(fecundity.summary$focal==focal)]),
                                                                  sd = log(fecundity.summary$st.dev.seed[which(fecundity.summary$focal==focal)]))), 
                                                            dim = 1))
                                 
  FinalFit <- stan(file = "code/DensityFunct_BH_Final.stan", 
                   data = DataVec,
                   init=  list.init,
                   warmup= 500,
                   iter = 1000, 
                   init_r = 2,
                   chains = 3,
                   control=list(max_treedepth=15),
                   seed= 165)
  
  
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





