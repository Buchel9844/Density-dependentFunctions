APM_PostFecunditydistribution  <- data.frame()
for(Code.focal in c("j","i")){

    # subset for one Code.focal species 
    SpDataFocal <- simulated.data[which(simulated.data$focal == Code.focal),]
    
    SpDataFocal$fecundity[is.na(SpDataFocal$fecundity)] <- 0
    
    
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
                                                                "seeds.j","fecundity")]
    
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
    
    #---- 2.2. Interaction (HOIs) matrix of plant species with COMPETITORS ----
    
    # creation of a matrix of S by S of the interaction jk in HOIs_ijk for plants
    matrix_HOIs_plant <- list()
    for (i in 1:N){
      matrix_i <- matrix(nrow=S,ncol=S)
      for (n in 1:S) {
        for (m in 1:S) {
          if (m <= n){
            matrix_i[n,m] = 0
          }
          else{
            matrix_i[n,m] = SpMatrix[i,n]* SpMatrix[i,m]
          } 
        }
      }
      matrix_i[is.na(matrix_i)] <- 0
      matrix_HOIs_plant[[i]] <-  matrix_i
    }
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #---- 3. BAYES FIT----
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##---- 3.1. Set up summary interactions df and parameters ---- 
    run_estimation <- 1
    DataVec <- c("N", "S",
                 "Fecundity", "SpMatrix",
                 "Intra","run_estimation")
    
    ##---- 3.2. Run  final fit ----
    # Now run a fianl fit of the model to assess parameter 
    print("Final Fit beginning")
    
    #install.packages("codetools")
    library("codetools")
    options(mc.cores = parallel::detectCores())
    FinalFit <- stan(file = "code/DensityFunct_APM_Final.stan", 
                    data = DataVec,
                    init=0,
                    warmup= 500,
                    iter = 1000, 
                    chains = 3)
    
    
    save(file= paste0("results/FinalFit_APM_",Code.focal,".rds"),
         FinalFit)

   load(file= paste0("results/FinalFit_APM_",Code.focal,".rds"))
    FinalPosteriors <- rstan::extract(FinalFit)
    
    print("Final Fit done")
    
    #---- 3.3. Final fit posterior check and behavior checks---- 
    
    ##### Diagnostic plots and post prediction 
    pdf(paste0("figures/FinalFit_APM_",Code.focal,".pdf"))
    # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
    source("code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
    # check the distribution of Rhats and effective sample sizes 
    ##### Posterior check
    # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
    hist(summary(FinalFit)$summary[,"Rhat"],
         main = paste("Finat Fit: Histogram of Rhat for",
                      Code.focal))
    hist(summary(FinalFit)$summary[,"n_eff"],
         main = paste("Finat Fit: Histogram of Neff for",
                      Code.focal))
    
    # plot the corresponding graphs
    stan_model_check(FinalFit,
                     param =c('lambdas','alpha_function_tilde'))
    
    # Next check the correlation among key model parameters and identify any
    pairs(FinalFit, pars = c("lambdas", 'alpha_function_tilde'))
    
    
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
                                      alpha_mean = -log(-Alphadistribution.i$mean),
                                      alpha_sd= -log(-Alphadistribution.i$sd),
                                      neighbours= "species i",focal = paste("species",Code.focal),
                                      density.function ="APM")
    
    
    Alphadistribution.j <- Alphadistribution.j %>% 
      group_by(species.j) %>% summarise_at("alpha.j",  list(mean = mean, sd = sd))
    
    Alphadistribution.j <- data.frame(abundance.neighbours = Alphadistribution.j$species.j,
                                      alpha_mean = -log(-Alphadistribution.j$mean),
                                      alpha_sd= -log(-Alphadistribution.j$sd),
                                      neighbours= "species j",focal = paste("species",Code.focal),
                                      density.function = "APM")
    
    
    Alphadistribution.neighbours <- bind_rows(Alphadistribution.neighbours,Alphadistribution.i,Alphadistribution.j)
    
    #---- 3.4. Extraction fecundity---
    
    APM_Fecunditydistribution.n <- FinalPosteriors %>% 
      as.data.frame() %>% 
      dplyr::select(contains("F_sim")) %>%
      gather( key="obervation",value="Fec")
    
    mu_apm <- FinalPosteriors[["F_sim"]]
    APM_Fecunditydistribution.n$obs <- rep(1:dim(mu_apm)[1],each=dim(mu_apm)[2])
    APM_Fecunditydistribution.n$iterations <- rep(1:dim(mu_apm)[2],times=dim(mu_apm)[1])
    
    APM_Fecunditydistribution.n$focal = Code.focal
    APM_Fecunditydistribution.n$alpha.function = "APM"
    
    APM_PostFecunditydistribution <- bind_rows(APM_PostFecunditydistribution,
                                               APM_Fecunditydistribution.n)
    
    
    
  }
write_csv(APM_PostFecunditydistribution, file = "results/APM_PostFecunditydistribution.csv")

