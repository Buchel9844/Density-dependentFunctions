# Function to create a simulated dataset 

  S = 2  # number of focal groups / species
  K = 2  # number of neighbour focals 
  pI = 0.1  # proportion of interactions which are NOT observed 

  
  # Neighbourhood abundance counts
  #-------------------------------
  # we assume we have the same number of observations for each focal
  S_obs <- c(100,100)
  K_Nmat <- NULL
  for(s in 1:S){ 
    # create a matrix of observations for each focal species
    foo <- matrix(data = 0, nrow = S_obs[s], ncol = K)
    for(j in 1:K) {
        for (i in 1:(S_obs[s]*1.5)){ 
          # randomly select an observation
          cellN <- round(runif(1, 1, S_obs[s]))
          # fill it with a 1 
          foo[cellN, j] <- foo[cellN, j] + 1  
          # and so on until neighbours are all accounted for
        }
    }
    if(identical(S_obs*1.5, colSums(foo)) == F) message('Error in abundance tallies!')
    K_Nmat <- rbind(K_Nmat, foo)
  }
  
  if (nrow(K_Nmat) != sum(S_obs)) message('Error in total number of observations!')
  colnames(K_Nmat) <- c("Spi","Spj") # species i invade species j
  
  # Parameters
  #-----------
  # log of species-specific intrinsic performance
  sim_a <- c(6,3)#order(runif(S, 2, 7)) # species i always has the higher one 
  
  # 'true' interaction strengths matrix as 
  #ii - ij 
  #ji - jj
  # intraspecific interactions negative and equal here
  # interspecific interaction of j on i is weaker than i on j
   sim_truealpha <- matrix(data = c(-0.1,-0.1,
                                    -0.3,-0.1),
                          nrow = S, ncol = K)
   
   # germination rate
   sim_g <- c(1,1)
   
   # seed survival rate
   sim_s <- c(1,1)
   
  # Observed seed set
  #------------------
  seeds <- rep(NA, length = sum(S_obs))
  counter <- 0
  for(s in 1:S) {
    for (i in (counter+1):(counter + S_obs[s])) {
      # multiply neighbour abundances by 'true' alpha values
      seeds[i] <- round(exp(sim_a[s] + sum(K_Nmat[i, ] * sim_g[s] * sim_truealpha[s, ])))
    }
    counter <- counter + S_obs[s]
  } 
  
  # Simulated dataset
  #------------------
  focal <- rep(c("Spi","Spj"), each=100)   # focal identifier
  fecundity <- as.numeric(seeds)
  simdata <- cbind(focal, fecundity, K_Nmat)
  simdata <- as.data.frame(simdata)








