# Recovery of parameters from simulation
library(cli)
library(broom)
library(colorspace)
library(truncnorm)
library(rstan)
library(loo) # Efficient approximate leave-one-out cross-validation
library(HDInterval)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(ggplotify) 
#library(odeintr)
#library(brms)
#library(bbmle)
#library(reshape2)
#library(ggthemes)
#library(data.table) # write csv.gz
#library(ppcor) # partial correlation value
#library(scales) # for pretty breaks with scale_axis
#library(grid)
#library(truncnorm)

home.dic <- "/home/lbuche/Eco_Bayesian/Complexity_caracoles/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"
#home.dic <-""

# ---- 1.0. Creates response surface data ----
source(paste0(home.dic,"code/stan/1.1.PopProjection_toolbox.R"))
#Fixed parameters
create.list.parameter=0
if(create.list.parameter ==1){
  params_torec.df <- NULL
  params_torec <- list()
  for(i in 1:100){
  params_torec[[i]] <- list(g = c(0.9,0.9),
                       s = c(0.9,0.9),
                       lambda = floor(abs(truncnorm::rtruncnorm(2, mean=20,sd=10,a=10,b=30))),
                       Nmax = matrix(ncol=2,nrow=2,
                                     floor(abs(truncnorm::rtruncnorm(4, mean=0,sd=3,a=0,b=10))),
                                     dimnames = list(c("i", "j"),
                                                     c("i", "j"))),
                       a_initial = matrix(ncol=2,nrow=2,
                                          truncnorm::rtruncnorm(4, mean=-0.05,sd=0.1,a=-0.5,b=0.5),
                                          dimnames = list(c("i", "j"),
                                                          c("i", "j"))),
                       c = matrix(ncol=2,nrow=2,
                                  truncnorm::rtruncnorm(4, mean=-0.05,sd=0.2,a=-0.5,b=0),
                                  dimnames = list(c("i", "j"),
                                                  c("i", "j"))),
                       a_slope = matrix(ncol=2,nrow=2,
                                        truncnorm::rtruncnorm(4, mean=-0.2,sd=0.2,a=-1,b=0),
                                        dimnames = list(c("i", "j"),
                                                        c("i", "j"))))
  
  params_torec.df.i <- params_torec[[i]] %>%
    as.data.frame() %>%
    rownames_to_column(var="focal")  %>%
    mutate(sim= i)
  names(params_torec.df.i) <- c("focal","g","s","lambda","N_opt_i","N_opt_j",
                              "alpha_init_i", "alpha_init_j",
                              "alpha_c_i","alpha_c_j",
                              "alpha_slope_i","alpha_slope_j","sim")
  params_torec.df <- bind_rows(params_torec.df,params_torec.df.i)
  }
  save(params_torec,
       file=paste0(home.dic,"results/stan/params_torec.rdata"))
  write.csv(params_torec.df,
        paste0(home.dic,"results/stan/params_torec.df.csv"))
}else{
  load(paste0(home.dic,"results/stan/params_torec.rdata"))
}
# percentage of obs equal zero
args <- commandArgs(trailingOnly = TRUE)
function.int.i <- as.numeric(args[1])

perc.zer <- c(0.2,0.4,0.5)
n.obs.levels <- c(50,150,300) #50
df.ite.param <- NULL

df.qtl.param <- NULL

#n.obs=100
#perc.zero.obs=1
#function.int=3
#it.bootstrap =89
for(it.bootstrap in 1:100){
  params_torec.i <- params_torec[[it.bootstrap]] 
  params_torec.df <- params_torec[[it.bootstrap]] %>%
    as.data.frame() %>%
    rownames_to_column(var="focal")
  names(params_torec.df) <- c("focal","g","s","lambda","N_opt_i","N_opt_j",
                                "alpha_init_i", "alpha_init_j",
                                "alpha_c_i","alpha_c_j",
                                "alpha_slope_i","alpha_slope_j")
for (n.obs in n.obs.levels){
  # initial density state
  for( perc.zero.obs in 1:3){
    # ---- 2.0. fit Bayesian model ----
    for (function.int in function.int.i){

        state.obs.i <- c(abs(round(rnorm(1000,5,10))))
        state.obs.j <- c(abs(round(rnorm(1000,5,10))))
        
        state.obs.i.n <- state.obs.i
        state.obs.j.n <- state.obs.j
        
        
        ts_param <-  as.data.frame(Ricker_solution_withalpha(gens=11,
                                                             state=c(state.obs.i.n[1],state.obs.j.n[1]),
                                                             pars = params_torec.i,
                                                             function.int=4,
                                                             add_external_factor="none")[1,])
        
        
        
        for(n in 2:1000){
          
          ts_param <- bind_rows(ts_param,as.data.frame(Ricker_solution_withalpha(gens=11,
                                                                                 state=c(state.obs.i.n[n],state.obs.j.n[n]),
                                                                                 pars =  params_torec.i,
                                                                                 function.int = 4,
                                                                                 add_external_factor="none")[1,]))
        }
        # remove a percentage of the observation equal to zero
        head(ts_param)
        #mean(ts_param$Fi)
        #mean(ts_param$Fj)
        
        vec.no.0 <- which(round(ts_param$Fi)>0 & round(ts_param$Fj)>0)
        if(length(vec.no.0) < n.obs){ next }
        
        fecundity.0 <- function(ts_param,n.obs,perc.zero.obs){
          repeat {
            # do something
            vec.no.0 <- which(round(ts_param$Fi)>0 & round(ts_param$Fj)>0)
            ts_param.n <- ts_param[vec.no.0[sample(1:length(vec.no.0),n.obs)],]
            F.zero.i <- sum(round(ts_param.n$Fi)==0)/n.obs
            F.zero.j <- sum(round(ts_param.n$Fj)==0)/n.obs
            # exit if the condition is met
            if (perc.zer[perc.zero.obs] >= F.zero.i &
                perc.zer[perc.zero.obs] >= F.zero.j) break
          }
          return(ts_param.n)
        }
        
        ts_param <- fecundity.0(ts_param,n.obs,perc.zero.obs)
        
        F.zero.i <- sum(round(ts_param$Fi)==0)/n.obs
        F.zero.j <- sum(round(ts_param$Fj)==0)/n.obs
        
        if (perc.zer[perc.zero.obs]-F.zero.i > 0){
          ts_param$Fi[sample(size=(perc.zer[perc.zero.obs]-F.zero.i )*n.obs,1:n.obs)] <- 0
        }
        
        if (perc.zer[perc.zero.obs]-F.zero.j > 0){
          ts_param$Fj[sample(size=(perc.zer[perc.zero.obs]-F.zero.j )*n.obs,1:n.obs)] <- 0
        }
        #plot(density( ts_param$Fj))
        #plot(density( ts_param$Fi))
        for (focal in c("i","j")){
          print(paste(n.obs,perc.zero.obs,function.int,it.bootstrap," for ",focal))
          
          function.vec <- c(0,0,0,0)
          function.vec[function.int] <- 1
          alphaFunct1 <- function.vec[1]
          alphaFunct2 <- function.vec[2]
          alphaFunct3 <- function.vec[3]
          alphaFunct4 <- function.vec[4]
          DataVec <- list(N=nrow(ts_param), 
                          n.obs=n.obs,
                          S=2,
                          lambda_prior = mean(ts_param[,paste0("F",focal)] + sd(ts_param[,paste0("F",focal)])),
                          Nmedian = c(round(unlist(c(ts_param[which.max(ts_param[,paste0("F",focal)]),
                                                                c("Stem_i","Stem_j")])))),# max fecundity
                          Fecundity=c(ceiling(ts_param[,paste0("F",focal)])),
                          SpMatrix =data.frame(Ni=ts_param$Stem_i,
                                               Nj=ts_param$Stem_j),
                          run_estimation=1,
                          alphaFunct1=alphaFunct1,
                          alphaFunct2=alphaFunct2,
                          alphaFunct3=alphaFunct3,
                          alphaFunct4=alphaFunct4
          )
          
          
          if(it.bootstrap == 1){
            
            assign(paste0("Parameters_",focal,"_",perc.zero.obs,"_",n.obs,"_",function.int),
                   list(DataVec = DataVec,
                        focal =focal,
                        params_torec =  params_torec.i,  
                        ts_param = ts_param
                   ))
            
            save(list =paste0("Parameters_",focal,"_",perc.zero.obs,"_",n.obs,"_",function.int),
                 file=paste0(home.dic,"results/stan/Parameters_",focal,"_",perc.zero.obs,"_",n.obs,"_",function.int,".RData"))
            
          }
          
          list.init <- function(...)list(lambda = array(as.numeric(DataVec$lambda_prior), 
                                                        dim = 1))
          
          rstan_options(auto_write = TRUE) 
          options(mc.cores = parallel::detectCores()) # to use the core at disposition 
          
          FinalFit <- rstan::stan(file = paste0(home.dic,"code/stan/CopyOfDensityFunct_Final.stan"), 
                           data = DataVec,
                           init=  list.init,
                           warmup= 500,
                           iter = 2000, 
                           init_r = 1,
                           chains = 3,
                           control=list(max_treedepth=15),
                           seed= 1644)
          
          FinalPosteriors <- rstan::extract(FinalFit)
          
          if(it.bootstrap == 1){
            print(paste(it.bootstrap ))
            save(file= paste0(project.dic ,"results/stan/FinalFit_",focal,"_",perc.zero.obs,"_",n.obs,"_",function.int,".rds"),
                 FinalFit)
            pdf(paste0(project.dic ,"figure/stan/FinalFit_",focal,"_",perc.zero.obs,"_",n.obs,"_",function.int,".pdf"))
            # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
            source(paste0(home.dic,"code/stan/stan_modelcheck_rem.R")) # call the functions to check diagnistic plots
            # check the distribution of Rhats and effective sample sizes 
            ##### Posterior check
            stan_post_pred_check(FinalPosteriors,"F_hat",
                                 DataVec$Fecundity)
            
            # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
            hist(summary(FinalFit)$summary[,"Rhat"],
                 main = paste("Finat Fit: Histogram of Rhat for",focal,"_",perc.zero.obs,"_",n.obs,
                              " and ",function.int))
            hist(summary(FinalFit)$summary[,"n_eff"],
                 main = paste("Finat Fit: Histogram of Neff for",focal,"_",perc.zero.obs,"_",n.obs,
                              " and ",function.int))
            
            # plot the corresponding graphs
            param.fit <- c('lambda',"N_opt",
                           'alpha_initial','alpha_slope','c',
                           "disp_dev")
            trace <- stan_trace(FinalFit, pars=param.fit,
                                inc_warmup = TRUE)
            print(trace)
            dens <- stan_dens(FinalFit, 
                              pars=param.fit)
            print(dens)
            splot <- stan_plot(FinalFit, 
                               pars=param.fit)
            print(splot)
            
            # Next check the correlation among key model parameters and identify any
            pairs(FinalFit, pars =param.fit)
            
            
            dev.off()
          }
          # ---- 2.1. Extract coefficients ----
          
          df.ite.param.n <- data.frame(focal =focal,
                                       lambda =  median(FinalPosteriors$lambda) - params_torec.df[which(params_torec.df$focal ==focal),
                                                                                                  "lambda"],
                                       N_opt_i =  median(FinalPosteriors$N_opt[,1])- params_torec.df[which(params_torec.df$focal ==focal),
                                                                                                     "N_opt_i"],
                                       N_opt_j =  median(FinalPosteriors$N_opt[,2])- params_torec.df[which(params_torec.df$focal ==focal),
                                                                                                     "N_opt_j"],
                                       alpha_slope_i = median(FinalPosteriors$alpha_slope[,1])- params_torec.df[which(params_torec.df$focal ==focal),
                                                                                                                "alpha_slope_i"],
                                       alpha_slope_j = median(FinalPosteriors$alpha_slope[,2])- params_torec.df[which(params_torec.df$focal ==focal),
                                                                                                                "alpha_slope_j"],
                                       alpha_init_i = median(FinalPosteriors$alpha_initial[,1])- params_torec.df[which(params_torec.df$focal ==focal),
                                                                                                                 "alpha_init_i"],
                                       alpha_init_j = median(FinalPosteriors$alpha_initial[,2])- params_torec.df[which(params_torec.df$focal ==focal),
                                                                                                                 "alpha_init_j"],
                                       alpha_c_i = median(FinalPosteriors$c[,1])- params_torec.df[which(params_torec.df$focal ==focal),
                                                                                                  "alpha_c_i"],
                                       alpha_c_j = median(FinalPosteriors$c[,2])- params_torec.df[which(params_torec.df$focal ==focal),
                                                                                                  "alpha_c_j"],
                                       perc.zero.obs=perc.zero.obs,
                                       n.obs=n.obs,function.int=function.int,
                                       it.bootstrap = it.bootstrap,
                                       mean.F = mean(DataVec$Fecundity),
                                       mean.N.i = mean(DataVec$SpMatrix$Ni),
                                       mean.N.j = mean(DataVec$SpMatrix$Nj))
          df.ite.param <- bind_rows(df.ite.param, df.ite.param.n)
          
          write.csv(df.ite.param,
                    paste0(home.dic,"results/stan/Parameters_median_",function.int,".csv"))
          
          mat.qtl.n <- matrix(c(quantile(FinalPosteriors$lambda,c(0.025,0.975)),
                                quantile(FinalPosteriors$N_opt[,1],c(0.025,0.975)),
                                quantile(FinalPosteriors$N_opt[,2],c(0.025,0.975)),
                                quantile(FinalPosteriors$alpha_slope[,1],c(0.025,0.975)),
                                quantile(FinalPosteriors$alpha_slope[,1],c(0.025,0.975)),
                                quantile(FinalPosteriors$alpha_initial[,1],c(0.025,0.975)),
                                quantile(FinalPosteriors$alpha_initial[,2],c(0.025,0.975)),
                                quantile(FinalPosteriors$c[,1],c(0.025,0.975)),
                                quantile(FinalPosteriors$c[,1],c(0.025,0.975))),
                              ncol=2, byrow = T)
          
          df.qtl.param.n <- data.frame(focal =focal,
                                       parameters.name = c("lambda",
                                                           "N_opt_i",
                                                           "N_opt_j","alpha_slope_i","alpha_slope_j",
                                                           "alpha_init_i","alpha_init_j",
                                                           "alpha_c_i","alpha_c_j"),
                                       Q2.5 =  mat.qtl.n[,1],
                                       Q97.5 =  mat.qtl.n[,2],
                                       perc.zero.obs=perc.zero.obs,
                                       n.obs=n.obs,function.int=function.int,
                                       it.bootstrap = it.bootstrap)
          df.qtl.param <- bind_rows(df.qtl.param, df.qtl.param.n)
          
          write.csv(df.qtl.param,
                    paste0(home.dic,"results/stan/Parameters_qtl_",function.int,".csv"))
        }
      }
    }
  }
}

write.csv(df.qtl.param,
          paste0(home.dic,"results/stan/Parameters_qtl_",function.int,".csv"))

write.csv(df.ite.param,
          paste0(home.dic,"results/stan/Parameters_median_",function.int,".csv"))

# ---- 3.0. Compare true and estimated median values ----

graph.output =0 
if(graph.output==1){
  Parameters_qtl <- read.csv(paste0(home.dic,"results/stan/Parameters_qtl_1.csv"))%>%
    bind_rows(read.csv(paste0(home.dic,"results/stan/Parameters_qtl_2.csv"))) %>%
    bind_rows(read.csv(paste0(home.dic,"results/stan/Parameters_qtl_3.csv"))) %>%
    bind_rows(read.csv(paste0(home.dic,"results/stan/Parameters_qtl_4.csv")))
  df.ite.param <- read.csv(paste0(home.dic,"results/stan/Parameters_median_1.csv")) %>%
    bind_rows(read.csv(paste0(home.dic,"results/stan/Parameters_median_2.csv"))) %>%
    bind_rows(read.csv(paste0(home.dic,"results/stan/Parameters_median_3.csv"))) %>%
    bind_rows(read.csv(paste0(home.dic,"results/stan/Parameters_median_4.csv")))
write.csv(Parameters_qtl,
          paste0(home.dic,"results/stan/Parameters_qtl.csv"))

write.csv(df.ite.param,
          paste0(home.dic,"results/stan/Parameters_median.csv"))

# percentage of obs equal zero
Parameters_qtl <- read.csv(paste0(home.dic,"results/stan/Parameters_qtl.csv")) 
df.ite.param <- read.csv(paste0(home.dic,"results/stan/Parameters_median.csv"))
perc.zer <- c(0.1,0.25,0.4)
n.obs.levels <- c(50,150,300)



param.fctlist <- list(param.fct.1 = c("lambda","alpha_init_i", "alpha_init_j"),
                      param.fct.2 = c("lambda","N_opt_i","N_opt_j",
                                      "alpha_init_i", "alpha_init_j",
                                      "alpha_slope_i","alpha_slope_j"),
                      param.fct.3 = c("lambda","N_opt_i","N_opt_j",
                                      "alpha_init_i", "alpha_init_j",
                                      "alpha_slope_i","alpha_slope_j",
                                      "alpha_c_i","alpha_c_j"),
                      param.fct.4 = c("lambda","N_opt_i","N_opt_j",
                                      "alpha_init_i", "alpha_init_j",
                                      "alpha_slope_i","alpha_slope_j",
                                      "alpha_c_i","alpha_c_j")
)

df.ite.param.n <- df.ite.param %>%
      #dplyr::filter(function.int == function.int ) %>%
      dplyr::select(focal,function.int,n.obs,perc.zero.obs, param.fctlist[[4]]) %>%
      gather(param.fctlist[[4]],
             key="parameter",value="median")  %>%
      mutate(parameter = case_when((parameter=="N_opt_i"|parameter=="N_opt_j")~ "N_opt",
                                   (parameter=="alpha_init_i"|parameter=="alpha_init_j")~ "alpha_init",
                                   (parameter=="alpha_slope_i"|parameter=="alpha_slope_j")~ "alpha_slope",
                                   (parameter=="alpha_c_i"|parameter=="alpha_c_j")~ "alpha_c",
                                   T~ parameter)) %>%
      filter((function.int == 1 & parameter %in% c("lambda","alpha_init"))|
               (function.int == 2 & parameter %in% c("lambda","alpha_init","alpha_slope","N_opt")) |
               (function.int == 3 & parameter %in% c("lambda","alpha_init","alpha_slope","alpha_c","N_opt")) |
               (function.int == 4 & parameter %in% c("lambda","alpha_init","alpha_slope","alpha_c","N_opt"))) %>%
  mutate(function.name = case_when(function.int==1 ~"Traditional constant",
                                   function.int==2 ~"Linear",
                                   function.int==3 ~"Exponential",
                                   function.int==4 ~"Sigmoid")) %>%
  mutate(function.name  = factor(function.name,
                                 levels=c("Traditional constant",
                                          "Linear","Exponential",
                                          "Sigmoid")))
my_cols <- c("#AA4499", "#DDCC77","#88CCEE", "#44AA99")   
hist.param <- df.ite.param.n %>%
  filter((parameter =="lambda" & median <25)| # for the xaxis to be symetric, does nto change the distributions
           (parameter =="N_opt" & median <10)|
           parameter %in% c("alpha_init","alpha_slope","alpha_c")) %>%
  mutate(parameter = factor(parameter,
                            levels=c("lambda","alpha_init","alpha_slope","N_opt","alpha_c"))) %>%
  ggplot(aes(x=median, # already removed the parameter to recover
             y=function.name,
             fill=as.factor(function.int))) +
  geom_density_ridges(quantile_lines = TRUE,
                      scale=1, 
                      quantiles = c(0.025,0.5, 0.975)) +
  scale_fill_manual(values=my_cols) +
  guides(fill="none") +
  facet_grid(.~parameter,scales = "free") +
  geom_vline(xintercept=0, color="red") +
  labs(y= "Function form used in the performance model\nto estimate parameters",
       x="Difference between the median of the estimated parameter and its corresponding set value") +
  theme_bw() +
  theme( legend.key.size = unit(1, 'cm'),
         legend.position = "bottom",
         strip.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(),
         strip.text = element_text(size=16),
         legend.text=element_text(size=12),
         legend.title=element_text(size=12),
         #axis.ticks.x=element_blank(),
         axis.text.x= element_text(size=14, angle=66, hjust=1),
         axis.text.y= element_text(size=14, angle=66, hjust=1),
         axis.title.x= element_text(size=16),
         axis.title.y= element_text(size=16),
         title=element_text(size=16),
         plot.margin = unit(c(1,1,1,1), "cm"))
hist.param
ggsave( hist.param,
        file=paste0(home.dic,"figure/stan/hist.param.fct.pdf"))

hist.param.n.obs <- df.ite.param.n %>%
      filter((parameter =="lambda" & median <20)| # for the xaxis to be symetric, does nto change the distributions
               (parameter =="N_opt" & median <10)|
               parameter %in% c("alpha_init","alpha_slope","alpha_c")) %>%
      mutate(parameter = factor(parameter,
                                levels=c("lambda","alpha_init","alpha_slope","N_opt","alpha_c"))) %>%
      ggplot(aes(x=median, # already removed the parameter to recover
                 y=function.name,
                 fill=as.factor(function.int))) +
      geom_density_ridges(quantile_lines = TRUE,
                          scale=1,
                          quantiles = c(0.025,0.5, 0.975)) +
  scale_fill_manual(values=my_cols) +
  guides(fill="none") +
  facet_grid(n.obs~parameter,scales = "free") +
      geom_vline(xintercept=0, color="red") +
      labs(y= "Function form used in the performance model\nto estimate parameters",
           x="Difference between the median of the estimated parameter and its corresponding set value") +
      theme_bw() +
      theme( legend.key.size = unit(1, 'cm'),
             legend.position = "bottom",
             strip.background = element_blank(),
             panel.grid.minor = element_blank(),
             panel.grid.major.x = element_blank(),
             strip.text = element_text(size=16),
             legend.text=element_text(size=12),
             legend.title=element_text(size=12),
             #axis.ticks.x=element_blank(),
             axis.text.x= element_text(size=12, angle=0, hjust=1),
             axis.text.y= element_text(size=12, angle=20, hjust=1),
             axis.title.x= element_text(size=16),
             axis.title.y= element_text(size=16),
             title=element_text(size=16))
    hist.param.n.obs
    
    ggsave( hist.param.n.obs,
            file=paste0(home.dic,"figure/stan/hist.param.n.obs.fct.pdf"))
    
  hist.param.n.perc.zero <- df.ite.param.n %>%
      filter((parameter =="lambda" & median <20)| # for the xaxis to be symetric, does nto change the distributions
               (parameter =="N_opt" & median <10)|
               parameter %in% c("alpha_init","alpha_slope","alpha_c")) %>%
      mutate(parameter = factor(parameter,
                                levels=c("lambda","alpha_init","alpha_slope","N_opt","alpha_c"))) %>%
      mutate(perc.zero.obs = case_when(perc.zero.obs ==1 ~ "10% of zero",
                                       perc.zero.obs ==2 ~ "25% of zero",
                                       perc.zero.obs ==3 ~ "40% of zero")) %>%
        ggplot(aes(x=median, # already removed the parameter to recover
                   y=function.name,
                   fill=as.factor(function.int))) +
      geom_density_ridges(quantile_lines = TRUE,
                          scale=1,
                          quantiles = c(0.025,0.5, 0.975)) +
    scale_fill_manual(values=my_cols) +
    guides(fill="none") +
    facet_grid(as.character(perc.zero.obs) ~ parameter,
                 scales = "free") +
      geom_vline(xintercept=0, color="red") +
      labs(y= "Function form used in the performance model\nto estimate parameters",
           x="Difference between the median of the estimated parameter and its corresponding set value") +
    theme_bw() +
    theme( legend.key.size = unit(1, 'cm'),
           legend.position = "bottom",
           strip.background = element_blank(),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_blank(),
           strip.text = element_text(size=16),
           legend.text=element_text(size=12),
           legend.title=element_text(size=12),
           #axis.ticks.x=element_blank(),
           axis.text.x= element_text(size=12, angle=66, hjust=1),
           axis.text.y= element_text(size=12, angle=20, hjust=1),
           axis.title.x= element_text(size=16),
           axis.title.y= element_text(size=16),
           title=element_text(size=16))
  hist.param.n.perc.zero
  
     ggsave( hist.param.n.perc.zero,
            file=paste0(home.dic,"figure/stan/hist.param.n.perc.zero.fct.pdf"))
    

#---- Pvalue
library(dplyr)

for (row.n in 1:nrow(Parameters_qtl)){
  focal.n <- Parameters_qtl$focal[row.n]
  param.n <- Parameters_qtl$parameters.name[row.n]
  low.Q <- Parameters_qtl$Q2.5[row.n]
  high.Q <- Parameters_qtl$Q97.5[row.n]
  
 params_torec.df.i <- params_torec[[Parameters_qtl$it.bootstrap[row.n]]]  %>%
    as.data.frame() %>%
    rownames_to_column(var="focal")
  names(params_torec.df.i) <- c("focal","g","s","lambda","N_opt_i","N_opt_j",
                              "alpha_init_i", "alpha_init_j",
                              "alpha_c_i","alpha_c_j",
                              "alpha_slope_i","alpha_slope_j")
  
  param.to.recover.n <-  params_torec.df.i[which(params_torec.df.i$focal ==focal.n),param.n]
  Parameters_qtl$param.to.recover[row.n] <- param.to.recover.n
  Parameters_qtl$count <-1
  if(low.Q < param.to.recover.n &
     high.Q > param.to.recover.n ){
    Parameters_qtl$bin.p.val[row.n] <- 1
  }else{Parameters_qtl$bin.p.val[row.n] <- 0
  }
  if(!param.n %in% param.fctlist[[Parameters_qtl$function.int[row.n]]]){
    Parameters_qtl$bin.p.val[row.n] <- "to.removed"
  }
}
     head(Parameters_qtl)
Parameters_qtl <-  Parameters_qtl %>%
  dplyr::filter(!bin.p.val =="to.removed")%>%
  mutate(bin.p.val = as.numeric(bin.p.val)) %>%
  mutate(parameter = case_when((parameters.name =="N_opt_i"|parameters.name =="N_opt_j")~ "N_opt",
                               (parameters.name =="alpha_init_i"|parameters.name =="alpha_init_j")~ "alpha_init",
                               (parameters.name =="alpha_slope_i"|parameters.name =="alpha_slope_j")~ "alpha_slope",
                               (parameters.name =="alpha_c_i"|parameters.name =="alpha_c_j")~ "alpha_c",
                               T~ parameters.name))  %>%
  mutate()

write.csv(Parameters_qtl,
          paste0(home.dic,"results/stan/Parameters_qtl.csv"))

#view(Parameters_qtl)
Parameters_qtl.count <- Parameters_qtl %>% 
  aggregate(count ~ parameter + function.int + 
              perc.zero.obs + n.obs, function(x){sum(x)}) %>%
  left_join(Parameters_qtl %>% 
              aggregate(count ~ parameter +function.int + 
                          n.obs, function(x){sum(x)})%>% rename("Count.n.obs.fct"=count)) %>%
  left_join(Parameters_qtl %>% 
              aggregate(count ~ parameter  +function.int + 
                          perc.zero.obs , function(x){sum(x)})%>% rename("Count.perc.fct"=count)) %>%
  left_join(Parameters_qtl %>% 
              aggregate(count ~ parameter + function.int + count,
                        function(x){sum(x)})%>% rename("Count.fct"=count)) 
head(Parameters_qtl.count)

Parameters_qtl.sum <- Parameters_qtl %>% 
  aggregate(bin.p.val ~ parameter + function.int + 
              perc.zero.obs + n.obs, function(x){sum(x)}) %>%
  left_join(Parameters_qtl %>% 
              aggregate(bin.p.val ~ parameter +function.int + 
                          n.obs, function(x){sum(x)})%>% rename("Sum.n.obs.fct"=bin.p.val)) %>%
  left_join(Parameters_qtl %>% 
              aggregate(bin.p.val ~ parameter  +function.int + 
                          perc.zero.obs , function(x){sum(x)})%>% rename("Sum.perc.fct"=bin.p.val)) %>%
  left_join(Parameters_qtl %>% 
              aggregate(bin.p.val ~ parameter + function.int ,
                        function(x){sum(x)})%>% rename("Sum.fct"=bin.p.val))  %>%
  left_join(Parameters_qtl.count) %>%
  mutate(Perc.fct = round(Sum.fct/Count.fct*100,digits=2),
         Perc.perc.fct = round(Sum.perc.fct/Count.perc.fct*100,digits=2),
         Perc.n.obs.fct = round(Sum.n.obs.fct/Count.n.obs.fct*100,digits=2),
         Perc.n.obs.perc.fct = round(bin.p.val/count*100,digits=2))
  
head(Parameters_qtl.sum)

Parameters_qtl.sum.fct <- Parameters_qtl.sum%>%
  dplyr::select(c("parameter","function.int","Perc.fct")) %>% 
  unique() %>%
  mutate(parameter = factor(parameter,
                            levels=c("lambda","alpha_init","alpha_slope","N_opt","alpha_c"))) %>%
  spread(parameter,Perc.fct) %>%
  mutate(function.int = case_when(function.int==1 ~"Traditional constant",
                                   function.int==2 ~"Linear",
                                   function.int==3 ~"Exponential",
                                   function.int==4 ~"Sigmoid")) 
Parameters_qtl.sum.fct
  write.csv(Parameters_qtl.sum.fct, file=paste0(home.dic,"results/stan/Parameters_qtl.sum.fct.csv"))
Parameters_qtl.sum.perc.fct <- Parameters_qtl.sum%>%
  dplyr::select(c("parameter","function.int","perc.zero.obs","Perc.perc.fct")) %>% 
  unique() %>%
  mutate(parameter = factor(parameter,
                            levels=c("lambda","alpha_init","alpha_slope","N_opt","alpha_c"))) %>%
  spread(parameter,Perc.perc.fct) %>%
  mutate(function.int = case_when(function.int==1 ~"Traditional constant",
                                   function.int==2 ~"Linear",
                                   function.int==3 ~"Exponential",
                                   function.int==4 ~"Sigmoid")) %>%
  mutate(perc.zero.obs = case_when(perc.zero.obs ==1 ~ "10% of zero",
                                   perc.zero.obs ==2 ~ "25% of zero",
                                   perc.zero.obs ==3 ~ "40% of zero")) 
Parameters_qtl.sum.perc.fct

write.csv(Parameters_qtl.sum.perc.fct,
          file=paste0(home.dic,"results/stan/Parameters_qtl.sum.perc.fct.csv"))


Parameters_qtl.sum.n.obs.fct <- Parameters_qtl.sum%>%
  dplyr::select(c("parameter","function.int","n.obs","Perc.n.obs.fct")) %>% 
  unique() %>%
  mutate(parameter = factor(parameter,
                            levels=c("lambda","alpha_init","alpha_slope","N_opt","alpha_c"))) %>%
  spread(parameter,Perc.n.obs.fct) %>%
  mutate(function.int= case_when(function.int==1 ~"Traditional constant",
                                   function.int==2 ~"Linear",
                                   function.int==3 ~"Exponential",
                                   function.int==4 ~"Sigmoid"))
write.csv(Parameters_qtl.sum.n.obs.fct,
          file=paste0(home.dic,"results/stan/Parameters_qtl.sum.n.obs.fct.csv"))

Parameters_qtl.sum.n.obs.fct

}
