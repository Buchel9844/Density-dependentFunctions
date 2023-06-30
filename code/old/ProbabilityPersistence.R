# Probability of persisting together 
#### with function 1 ####
# distribution is uniform, so persisting probability is drwn from a uniform  
#### Extract lambdas ####
nspec <- 2 #number of species
lambda.df <- data.frame(alpha.function =rep(c("function_1","function_2",
                                              "function_3","function_4"), times=nspec),
                        focal = rep(c("species i","species j"),each=4))
for(alpha.function in c("function_1","function_2",
                        "function_3","function_4")){
  for(focal in c("i","j")){
load(paste0("results/FinalFit_",focal,"_",alpha.function,".rds"))
FinalPosteriors <- rstan::extract(FinalFit)

lambda.df[which(lambda.df$focal== paste0("species ",focal) &
                  lambda.df$alpha.function ==alpha.function),"mean.lambda"] <- mean(FinalPosteriors$lambdas)
lambda.df[which(lambda.df$focal== paste0("species ",focal) &
                  lambda.df$alpha.function ==alpha.function),"sd.lambda"] <- sd(FinalPosteriors$lambdas)

remove(FinalFit)
remove(FinalPosteriors)
  }
}
#check that the lambda follows a normal distribution:
ggplot(as.data.frame(FinalPosteriors$lambdas)) + geom_density(aes(x=V1))

#### Extract alpha matrix ####

df.proba <- NULL

for(i in 1:100){
set.seed(i) # to insure reproductible results
for(alpha.function in c("function_1","function_2",
                        "function_3","function_4")){
  alpha.matrix <- matrix(nrow=2,ncol=2,
                         dimnames = list(c("species i","species j"),
                                         c("species i","species j")))
  r.vector <- data.frame(focal = c("species i","species j"),
                         lambda=NA)
for(focal in c("species i","species j")){
  
  for(comp in c("species i","species j")){
    alpha.df <- Alphadistribution.neighbours[which(Alphadistribution.neighbours$density.function == alpha.function &
                                                   Alphadistribution.neighbours$focal == focal &
                                                   Alphadistribution.neighbours$neighbours == comp),]
    alpha.matrix[focal,comp] <- rnorm(1,mean = mean(alpha.df$alpha_mean),
                          sd = sum(alpha.df$alpha_sd))
  }
  r.vector[which(r.vector$focal==focal),"lambda"] <- c(rnorm(1,mean = lambda.df[which(lambda.df$focal== focal &
                                                      lambda.df$alpha.function ==alpha.function),"mean.lambda"],
                             sd = lambda.df[which(lambda.df$focal== focal &
                                                    lambda.df$alpha.function ==alpha.function),"sd.lambda"])
                       )# order i,j
  
   }

  df.proba.n <- data.frame(time = i,
                           alpha.function = alpha.function,
                           alpha_ii = c(alpha.matrix)[1],
                              alpha_ji = c(alpha.matrix)[2],
                              alpha_ij = c(alpha.matrix)[3],
                              alpha_jj = c(alpha.matrix)[4],
                              r.lambda.i =  r.vector[which(r.vector$focal=="species i"),
                                                     "lambda"],
                           r.lambda.j =  r.vector[which(r.vector$focal=="species j"),
                                                  "lambda"],
                           ND = Omega(alpha.matrix),
                           centroid_i = r_centroid(alpha.matrix)[1,1],
                           centroid_i = r_centroid(alpha.matrix)[2,1],
                           FD = theta(alpha.matrix,r.vector$lambda),
                           proba = test_feasibility(alpha.matrix,r.vector$lambda),
                           overlap = compute_overlap(alpha.matrix,100)[["overlap"]])
  df.proba <- bind_rows(df.proba, df.proba.n)
  }
}

df.proba.recap <- aggregate(proba~ alpha.function,
                            df.proba,
                            function(x) sum(x))
