# This functions spits out diagnostics for convergence
stan_diagnostic <- function(fit,
                            ...) {
  
  # Print diagnostics
  check_hmc_diagnostics(fit)
  
  check_treedepth(fit)
  # if fit saturates the threshold, rerun with larger maximum tree depth and recheck saturation
  check_energy(fit) # check E-BFMI
  check_divergences(fit) # check divergence (validity concern)
  
  fit_summary <- summary(fit)$summary
  # get the range of rhat values
  print(c('rhat range: ', range(na.omit(fit_summary[ , 'Rhat']))))
  # and n_eff
  print(c('n_eff range: ', range(na.omit(fit_summary[ , 'n_eff']))))
  
  # print some diagnostic plots 
  x <- stan_rhat(fit)  # distribution of Rhat
  print(x)
  x <- stan_ess(fit)   # ratio of effective sample size to total sample size
  print(x)
  x <- stan_mcse(fit)  # ratio of Monte Carlo standard error to posterior standard deviation
  print(x)
  
  # FOLLOWING ONLY WORK ON MULTIPLE CHAINS
  # #Histograms of lp__ and accept_stat, and their joint distribution
  stan_diag(fit, 'sample')
  # # Violin plots showing the distributions of lp__ and accept_stat at each of the sampled step sizes (one per chain)
  stan_diag(fit, 'stepsize')
  # # Histogram of treedepth and violin plots showing the distributions of lp__ and accept_stat for each value of treedepth
  stan_diag(fit, 'treedepth')
  # # Violin plots showing the distributions of lp__ and accept_stat for iterations that encountered divergent transitions (divergent=1) and those that did not (divergent=0)
  stan_diag(fit, 'divergence')

  
  # NB: for further diagnostics, I can explore with
  # - stan_par(fit, 'specific parameter')
  # - stan_ac(fit, 'param') for auto-correlation
}



# This function checks the validity of the RE model

stan_model_check <- function(fit,
                             params,
                             ...) {
  
  fit_summary <- summary(fit)$summary
  
  # remove 'mu' from the parameter vector and do it after because it is too long
  params <- params[!params %in% c('mu', 'mu2', 'log_lik_rim', 'log_lik_nddm')]
  
  sapply(params, function(x) { 
    
    # this is just to length the plots so they can show all parameters
    N <- length(grep(x, rownames(fit_summary)))
    # some exceptions: 
    if (x == 'beta_i0') {N <- 20}
    if (x == 'beta_ij') {N <- 800}
    if (x == 'lp__') {N <- 10}
    
    # save traceplot of the posterior draws
    tplot <- traceplot(fit, pars = x, inc_warmup = F, ncol = 5)
  
    print(tplot)

    
    # plot posterior uncertainty intervals
    postint <- plot(fit, pars = x, show_density = T)
    print(postint)
    
    })
  
  
}


# This function plots the posterior predicted seed number vs the actual data

stan_post_pred_check_nbin <- function(post.draws,
                                 var_name,
                                 stan.data,
                                 file.name,
                                 limx,
                                 ...) {
  
  # currently using the loo package, can switch to rethinking
  
  # phi is the overdispersion parameter for the neg binom model
  # mu is the mean for predicted seed number 
  
  # extract mu and phi
  mu <- post.draws[[var_name]] # matrix with nrow = draws and ncol = observations
  disp_dev <- post.draws$disp_dev
  phi <- (disp_dev^2)^(-1)
  
  # generating posterior predictions
  seed_pred <- matrix(nrow = dim(mu)[1], ncol = dim(mu)[2])
  for (i in 1:dim(mu)[1]) {     # for each posterior draw
    for (j in 1:dim(mu)[2]) {    # for each observation 
      # draw from the predicted distribution
      seed_pred[i, j] <- rnbinom(1, mu = mu[i, j], size = phi[i])
      if(is.na(seed_pred[i, j])){
        seed_pred[i, j] <- 0
      }
      if(seed_pred[i, j] > 10000){
      seed_pred[i, j] <- 10000
      }
    }
  }
  
  seed_pred_table <- gather(as_tibble(seed_pred),key = "key",
                            value = "Fec")
  seed_pred_table$obs <- rep(1:dim(mu)[1],each=dim(mu)[2])
  seed_pred_table$iterations <- rep(1:dim(mu)[2],times=dim(mu)[1])
  # get maximum density for plot limits
  max.density <- max(c(apply(seed_pred, 1, function(x) {max(density(x,na.rm=T)$y)}), 
                       max(density(stan.data)$y)))
  # dev.new(noRStudioGD = T)
  # start a plot with the first draw 
  ppc.plot <- plot(density(seed_pred[1, ]), ylim = c(0, max.density), col = 'darkgrey',
                   xlim =c(0,limx),
                   ylab = 'Seed density',
                   main = 'Post. pred. check',
                   sub = '(grey = predicted, black = observed)') 
  for (i in 2:dim(seed_pred)[1]) {
    # add a line for each draw
    ppc.plot <- lines(density(seed_pred[i, ]), col = 'darkgrey')
  }
  # add the actual data
  ppc.plot <- lines(density(stan.data), col = 'black', lwd = 2)  
  print(ppc.plot)
  
  write_csv(seed_pred_table, file = file.name)
}

stan_post_pred_check_norm <- function(post.draws,
                                 var_name = 'mu',
                                 stan.data,
                                 ...) {
  #post.draws =Inv.ModelfitPosteriors
  #var_name = "seedtotal"
  #stan.data=DataVec$obs_spI
  # currently using the loo package, can switch to rethinking
  
  # phi is the overdispersion parameter for the neg binom model
  # mu is the mean for predicted seed number 
  
  # extract mu and phi
  mu <- post.draws[[var_name]]*post.draws$g_init[1] # matrix with nrow = draws and ncol = observations
  disp_dev <- post.draws$disp_dev
  
  # generating posterior predictions
  seed_pred <- matrix(nrow = dim(mu)[1], ncol = dim(mu)[2])
  for (i in 1:dim(mu)[1]) {     # for each posterior draw
    for (j in 1:dim(mu)[2]) {    # for each observation 
      # draw from the predicted distribution
      seed_pred[i, j] <- rlnorm(1, mean= mu[i, j], sd = 1/disp_dev[i,]^2)
    }
  }

  # get maximum density for plot limits
  max.density <- max(c(apply(seed_pred, 1, function(x) {max(density(x)$y)}),
                       max(density(stan.data)$y)))
   max.density.x <- ceiling(max(density(stan.data)$x) + max(density(stan.data)$x)/10)
  
  # dev.new(noRStudioGD = T)
  # start a plot with the first draw 
  ppc.plot <- plot(density(seed_pred[1, ]), 
                   #ylim = c(0, max.density),
                   xlim=c(0, max.density.x),
                   ylab = 'Seed density',
                   xlab="") 
  for (i in 2:dim(seed_pred)[1]) {
    # add a line for each draw
    ppc.plot <- lines(density(seed_pred[i, ]))
  }
  # add the actual data
  ppc.plot <- lines(density(stan.data), col="red",lwd = 3)  
  #ppc.plot <- text(labels = value.se, 
    #               y= 0.04,#max.density - max.density/10,
    #               x=50, #max(stan.data)- max(stan.data)/10,
     #              cex=1, pos=3,col="black")
  print(ppc.plot)
}


stan_post_pred_check_pois <- function(post.draws,
                                      var_name,
                                      stan.data,
                                      ...) {
  
  # currently using the loo package, can switch to rethinking
  
  # phi is the overdispersion parameter for the neg binom model
  # mu is the mean for predicted seed number 
  
  # extract mu and phi
  mu <- post.draws[[var_name]]*post.draws$g_init[1] # matrix with nrow = draws and ncol = observations
  
  # generating posterior predictions
  seed_pred <- matrix(nrow = dim(mu)[1], ncol = dim(mu)[2])
  for (i in 1:dim(mu)[1]) {     # for each posterior draw
    for (j in 1:dim(mu)[2]) {    # for each observation 
      # draw from the predicted distribution
      seed_pred[i, j] <- rpois(1, lambda= mu[i, j])  
    }
  }
  
    # get maximum density for plot limits
  max.density <- max(c(apply(seed_pred, 1, function(x) {max(density(x,na.rm=T)$y)}), 
                       max(density(stan.data)$y)))
  max.density.x <- max(density(stan.data)$x)
  
  # dev.new(noRStudioGD = T)
  # start a plot with the first draw 
  ppc.plot <- plot(density(seed_pred[1, ],na.rm=T), 
                   ylim = c(0, max.density), 
                   xlim=c(0, max.density.x),
                   col = 'darkgrey',
                   ylab = 'Seed density',
                   main = 'Post. pred. check',
                   sub = '(grey = predicted, black = observed)') 
  for (i in 2:dim(seed_pred)[1]) {
    # add a line for each draw
    ppc.plot <- lines(density(seed_pred[i, ],na.rm=T), col = 'darkgrey')
  }
  # add the actual data
  ppc.plot <- lines(density(stan.data), col = 'red', lwd = 2)  
  print(ppc.plot)
  
 # write_csv(seed_pred_table, file = file.name)
}


stan_post_pred_check_pois_glow_ghigh <- function(post.draws,
                                      var_name,
                                      stan.data,
                                      gvec,
                                      numb_plot,
                                      numb_year,
                                      ...) {
  
  # currently using the loo package, can switch to rethinking
  
  # phi is the overdispersion parameter for the neg binom model
  # mu is the mean for predicted seed number 

  # extract mu and phi
  mu <- post.draws[[var_name]]* rep(gvec[1:numb_year-1], each=numb_plot) # matrix with nrow = draws and ncol = observations
  
  # generating posterior predictions
  seed_pred <- matrix(nrow = dim(mu)[1], ncol = dim(mu)[2])
  for (i in 1:dim(mu)[1]) {     # for each posterior draw
    for (j in 1:dim(mu)[2]) {    # for each observation 
      # draw from the predicted distribution
      seed_pred[i, j] <- rpois(1, lambda= mu[i, j])  
    }
  }
  
  # get maximum density for plot limits
  max.density <- max(c(apply(seed_pred, 1, function(x) {max(density(x,na.rm=T)$y)}), 
                       max(density(stan.data)$y)))
  max.density.x <- max(density(stan.data)$x)
  
  # dev.new(noRStudioGD = T)
  # start a plot with the first draw 
  ppc.plot <- plot(density(seed_pred[1, ],na.rm=T), 
                   ylim = c(0, max.density), 
                   xlim=c(0, max.density.x),
                   col = 'darkgrey',
                   ylab = 'Seed density',
                   main = 'Post. pred. check',
                   sub = '(grey = predicted, black = observed)') 
  for (i in 2:dim(seed_pred)[1]) {
    # add a line for each draw
    ppc.plot <- lines(density(seed_pred[i, ],na.rm=T), col = 'darkgrey')
  }
  # add the actual data
  ppc.plot <- lines(density(stan.data), col = 'red', lwd = 2)  
  print(ppc.plot)
  
  # write_csv(seed_pred_table, file = file.name)
}





