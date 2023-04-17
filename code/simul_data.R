# Function to create a simulated dataset 
simul_data <- function(
    # number of species has to be a multiple of 10
  S = 2,   # number of focal groups / species
  K = 2,   # number of neighbour focals 
  pI = 0.1,   # proportion of interactions which are NOT observed 
  ...
) {
  # Neighbourhood abundance counts
  #-------------------------------
  # we assume we have a different number of observations for each focal
  S_obs <- rep(100,times=S)
  
  # we assume that S = K, and that the number of observations of K is relative to their abundance
  approx_K_obs <- round(jitter(10*(S_obs), amount = 30))
  summedSK <- round(jitter(sapply(S_obs, '*', approx_K_obs)/1000, 20))
  
  # set up an empty matrix with columns for each neighbours
    # create a matrix of observations for each focal species
  K_Nmat <- matrix(ncol = K)
  for(s in 1:S){ 
    # create a matrix of observations for each focal species
    foo <- matrix(data = 0, nrow = S_obs[s], ncol = K)
    for(j in 1:K) {
        for (i in 1:500){ 
          # randomly select an observation
          cellN <- round(runif(1, 1, S_obs))
          # fill it with a 1 
          foo[cellN, j] <-   foo[cellN, j] + 1  
          # and so on until neighbours are all accounted for
        }
    }
    K_Nmat <- rbind(K_Nmat, foo)
  }
  K_Nmat <- K_Nmat[-1,]
  colSums(K_Nmat)
  
  if (nrow(K_Nmat) != sum(S_obs)) message('Error in total number of observations!')
  colnames(K_Nmat) <- c("i","j")
  Spmax <- max(K_Nmat,na.rm=T)
  if(!Spmax < 15){print("start again spmax too low")}
  
  #K_Nmat <-round((K_Nmat/max(K_Nmat))*100) #scale all the abundance between 0 and 100
  
  # Parameters for within guilt pairwise interaction
  #-------------------------------------------------
  # log of species-specific intrinsic performance
  sim_lambda <- runif(S, 5, 10)  
  
  # 'true' interaction strengths based on random draw
  #a <- rnorm( 8000, 0, .3 )
  #b <- rnorm( 2000, -0.7, .3 )
  #c <- rnorm( 1000,  0.3, .2 )
  #customdist <- c( a, b, c ) ## make a weird dist with Kurtosis and Skew
  #sim_truealpha <- matrix(data = sample(customdist, S*K, replace = F),
  #                       nrow = S, ncol = K)
  
  # 'true' interaction strengths based choosen values

  sim_alphaspecific <- matrix(nrow=S,ncol=K,
                              sample(seq(-1,1,by=0.01),S*K, replace = TRUE))
  
  diag(sim_alphaspecific) <- - abs(diag(sim_alphaspecific))
  
  # Observed seed set
  #------------------
  seeds <- rep(NA, length = sum(S_obs))
  counter <- 1
  
  for(s in 1:S) {
    
    for (i in (counter):(counter + S_obs[s] - 1)) {
      # print(i)
      # multiply neighbour abundances by 'true' alpha values

        seeds[i] <- round(exp(sim_lambda[s] + sum(K_Nmat[i, ] * sim_alphaspecific[s, ])))
      }
    counter <- counter + S_obs[s]
  }
  
  # seeds <- round(log(seeds)) # to include facilitation we took the exponential of the above expression. 
  # To have the "observed" fecundity we need to take its natural logarithm
  # Simulated dataset
  #------------------
  focal <- do.call('c', mapply(rep, c("i","j"), S_obs, SIMPLIFY = F))    # focal identifier

  colnames(sim_alphaspecific) <- paste0('alpha', c("i","j"))
  row.names(sim_alphaspecific) <- paste0('focal', c("i","j"))

    simdata <- cbind(focal, seeds, K_Nmat)
    simdata <- as.data.frame(simdata)
    
    return(list(simdata= simdata, sim_lambda = sim_lambda,
                sim_alpha_specific=sim_alphaspecific))
}


