# Recovery of parameters from simulation
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
library(odeintr)
library(brms)
library(bbmle)
library(reshape2)
library(ggthemes)
library(data.table) # write csv.gz
library(ppcor) # partial correlation value
library(scales) # for pretty breaks with scale_axis
library(grid)

home.dic <- "/home/lbuche/Eco_Bayesian/Complexity_caracoles/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"


# ---- 1.0. Creates response surface data ----
source(paste0(home.dic,"code/stan/1.1.PopProjection_toolbox.R"))
#Fixed parameters

params_torec <- list(g = c(0.9,0.9),
                     s = c(0.9,0.9),
                     lambda = c(4,4),
                     Nmax = matrix(ncol=2,nrow=2,
                                   c(0,1,3,0),
                                   dimnames = list(c("i", "j"),
                                                   c("i", "j"))),
                     a_initial = matrix(ncol=2,nrow=2,
                                        c(-0.10,-0.05,0.05,-0.10),
                                        dimnames = list(c("i", "j"),
                                                        c("i", "j"))),
                     a_slope = matrix(ncol=2,nrow=2,
                                      c(-0.1,-0.05,-0.05,-0.1),
                                      dimnames = list(c("i", "j"),
                                                      c("i", "j"))),
                     c = matrix(ncol=2,nrow=2,
                                c(-0.3,-0.25,-0.3,-0.35),
                                dimnames = list(c("i", "j"),
                                                c("i", "j"))))

params_torec.df <- params_torec %>%
  as.data.frame() %>%
  rownames_to_column(var="focal") 
names(params_torec.df) <- c("focal","g","s","lambda","N_opt_i","N_opt_j",
                            "alpha_init_i", "alpha_init_j",
                            "alpha_slope_i","alpha_slope_j",
                            "alpha_c_i","alpha_c_j")

# percentage of obs equal zero
perc.zer <- c(0.2,0.4,0.5)
n.obs.levels <- c(50,100,200,500)
df.ite.param <- read.csv(paste0(home.dic,
                                "results/stan/Parameters_median.csv"))

df.qtl.param <- read.csv(paste0(home.dic,
                                "results/stan/Parameters_qtl.csv"))

#n.obs=50
#perc.zero.obs=1
#function.int=3
#it.bootstrap =89
for (n.obs in 500){
  # initial density state
  for( perc.zero.obs in 1:3){
    # ---- 2.0. fit Bayesian model ----
    for (function.int in c(1:4)){
      for(it.bootstrap in 1:100){
        state.obs.i <- c(abs(round(rnorm(1000,5,10))))
        state.obs.j <- c(abs(round(rnorm(1000,5,10))))
        
        state.obs.i.n <- state.obs.i
        state.obs.j.n <- state.obs.j
        
        
        ts_param <-  as.data.frame(Ricker_solution_withalpha(gens=2,
                                                             state=c(state.obs.i.n[1],state.obs.j.n[1]),
                                                             pars = params_torec,
                                                             function.int=4,
                                                             add_external_factor="none")[1,])
        
        
        
        for(n in 2:1000){
          
          ts_param <- bind_rows(ts_param,as.data.frame(Ricker_solution_withalpha(gens=2,
                                                                                 state=c(state.obs.i.n[n],state.obs.j.n[n]),
                                                                                 pars = params_torec,
                                                                                 function.int = 4,
                                                                                 add_external_factor="none")[1,]))
        }
        # remove a percentage of the observation equal to zero
        #mean(ts_param$Fi)
        #mean(ts_param$Fj)
        
        
        fecundity.0 <- function(ts_param,n.obs,perc.zero.obs){
          repeat {
            # do something
            ts_param.n <- ts_param[sample(1:1000,n.obs),]
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
                          Nmedian = c(unlist(c(ts_param[which.max(ts_param[,paste0("F",focal)]),
                                                        c("Ni","Nj")]))),# max fecundity
                          Fecundity=c(round(ts_param[,paste0("F",focal)])),
                          SpMatrix =data.frame(Ni=ts_param$Ni*params_torec$g[1],
                                               Nj=ts_param$Nj*params_torec$g[1]),
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
                        params_torec = params_torec,  
                        ts_param = ts_param
                   ))
            
            save(list =paste0("Parameters_",focal,"_",perc.zero.obs,"_",n.obs,"_",function.int),
                 file=paste0(home.dic,"results/stan/Parameters_",focal,"_",perc.zero.obs,"_",n.obs,"_",function.int,".RData"))
            
          }
          
          list.init <- function(...)list(lambda = array(as.numeric(round(max(log(DataVec$Fecundity))), 
                                                                   dim = 1))#,
                                         #N_opt = array(as.numeric(DataVec$Nmedian), 
                                         #                          dim = 2)
          )
          
          
          FinalFit <- stan(file = paste0(home.dic,"code/stan/DensityFunct_Final.stan"), 
                           data = DataVec,
                           init=  list.init,
                           warmup= 500,
                           iter = 1000, 
                           init_r = 1,
                           chains = 4,
                           control=list(max_treedepth=15),
                           seed= 1644)
          FinalPosteriors <- rstan::extract(FinalFit)
          
          if(it.bootstrap == 1){
            print(paste(it.bootstrap ))
            save(file= paste0(home.dic,"results/stan/FinalFit_",focal,"_",perc.zero.obs,"_",n.obs,"_",function.int,".rds"),
                 FinalFit)
            pdf(paste0(home.dic,"figure/stan/FinalFit_",focal,"_",perc.zero.obs,"_",n.obs,"_",function.int,".pdf"))
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
            #pairs(FinalFit, pars = c("lambdas",'alpha_initial','alpha_slope','c'))
            
            
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
                    paste0(home.dic,"results/stan/Parameters_median.csv"))
          
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
                    paste0(home.dic,"results/stan/Parameters_qtl.csv"))
        }
      }
    }
  }
}

write.csv(df.qtl.param,
          paste0(home.dic,"results/stan/Parameters_qtl.csv"))

write.csv(df.ite.param,
          paste0(home.dic,"results/stan/Parameters_median.csv"))

# ---- 3.0. Compare true and estimated median values ----

params_torec.df <- params_torec %>%
  as.data.frame() %>%
  rownames_to_column(var="focal") 
names(params_torec.df) <- c("focal","g","s","lambda","N_opt_i","N_opt_j",
                            "alpha_init_i", "alpha_init_j",
                            "alpha_slope_i","alpha_slope_j",
                            "alpha_c_i","alpha_c_j")


Parameters_qtl <- read.csv(paste0(home.dic,"results/stan/Parameters_qtl.csv"))
df.ite.param <- read.csv(paste0(home.dic,"results/stan/Parameters_median.csv"))

# percentage of obs equal zero
perc.zer <- c(0.1,0.25,0.4)
n.obs.levels <- c(50,100,200,500)


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

list.plot <- list()
for (function.int in c(1:4)){
  for( perc.zero.obs in 1){
    
    df.ite.param.n <- df.ite.param %>%
      dplyr::filter(perc.zero.obs==perc.zero.obs &
                      function.int == function.int) %>%
      dplyr::select(focal,n.obs, param.fctlist[[function.int]]) %>%
      gather(param.fctlist[[function.int]],
             key="parameter",value="median") 
    
    params_torec.df.n <- params_torec.df %>%
      dplyr::select(focal, param.fctlist[[function.int]])%>%
      gather(param.fctlist[[function.int]],
             key="parameter",value="true.value") %>%
      mutate(focal.int = case_when(focal=="i"~1,
                                   focal =="j"~2))
    
    
    hist.param.n <- df.ite.param.n %>%
      ggplot(aes(x=median,
                 y=focal)) +
      geom_density_ridges(quantile_lines = TRUE,
                          scale=1,
                          quantiles = c(0.025,0.5, 0.975)) +
      facet_wrap(n.obs~parameter,scales = "free_x",nrow=4) +
      geom_point(data=params_torec.df.n,
                 aes(x=true.value,
                     y=as.factor(focal))) +
      geom_segment(data = params_torec.df.n, 
                   aes(x = true.value, 
                       xend = true.value, 
                       y = focal.int,
                       yend = focal.int +0.95),
                   color = "red", inherit.aes = F,
                   size=1)
    list.plot[[function.int]] <- hist.param.n
    
    ggsave( hist.param.n,
            file=paste0(home.dic,"figure/stan/hist.param.n.fct.",
                        function.int,".pdf"))
    
  }
}
list.plot[[1]]
list.plot[[2]]
list.plot[[3]]
list.plot[[4]]

#---- Pvalue
library(dplyr)


Parameters_qtl <- read.csv(paste0(home.dic,"results/stan/Parameters_qtl.csv"))

for (row.n in 1:nrow(Parameters_qtl)){
  focal.n <- Parameters_qtl$focal[row.n]
  param.n <- Parameters_qtl$parameters.name[row.n]
  low.Q <- Parameters_qtl$Q2.5[row.n]
  high.Q <- Parameters_qtl$Q97.5[row.n]
  param.to.recover.n <-  params_torec.df[which(params_torec.df$focal ==focal.n),param.n]
  Parameters_qtl$param.to.recover[row.n] <- param.to.recover.n
  if(low.Q < param.to.recover.n &
     high.Q > param.to.recover.n ){
    Parameters_qtl$bin.p.val[row.n] <- 1
  }else{Parameters_qtl$bin.p.val[row.n] <- 0
  }
  if(!param.n %in% param.fctlist[[Parameters_qtl$function.int[row.n]]]){
    Parameters_qtl$bin.p.val[row.n] <- "to.removed"
  }
}
Parameters_qtl <-  Parameters_qtl %>%
  dplyr::filter(!bin.p.val =="to.removed")%>%
  mutate(bin.p.val = as.numeric(bin.p.val))

write.csv(Parameters_qtl,
          paste0(home.dic,"results/stan/Parameters_qtl.csv"))

view(Parameters_qtl)
Parameters_qtl %>%
  mutate(parameters.name = factor(parameters.name,
                                  levels=c("lambda",
                                           "alpha_init_i", "alpha_init_j",
                                           "alpha_slope_i","alpha_slope_j",
                                           "alpha_c_i","alpha_c_j",
                                           "N_opt_i","N_opt_j"))) %>%
  ggplot(aes(group=function.int,x=bin.p.val)) +
  stat_count(aes(group=function.int,fill=as.factor(function.int))) +
  facet_grid(parameters.name~focal) +
  theme_bw()

Parameters_qtl.sum <- Parameters_qtl %>% 
  aggregate(bin.p.val ~ parameters.name + focal +function.int + 
              perc.zero.obs + n.obs, function(x) sum(x)/100 ) %>%
  left_join(Parameters_qtl %>% 
              aggregate(bin.p.val ~ parameters.name + focal +function.int + 
                          n.obs, function(x) sum(x)/300 )%>% rename("Sum.n.obs.fct"=bin.p.val)) %>%
  left_join(Parameters_qtl %>% 
              aggregate(bin.p.val ~ parameters.name + focal +function.int,
                        function(x) sum(x)/1200 )%>% rename("Sum.fct"=bin.p.val))

view(Parameters_qtl.sum)

ProptoRecoverParam <- Parameters_qtl.sum %>%
  mutate(parameters.name = factor(parameters.name,
                                  levels=c("lambda",
                                           "alpha_init_i", "alpha_init_j",
                                           "alpha_slope_i","alpha_slope_j",
                                           "alpha_c_i","alpha_c_j",
                                           "N_opt_i","N_opt_j"))) %>%
  ggplot(aes(x=as.factor(function.int),
             y=bin.p.val)) +
  geom_boxplot(aes(group=function.int,fill=as.factor(function.int))) +
  facet_wrap( #.~parameters.name,
    #n.obs~parameters.name,nrow=4, #perc.zero.obs
    perc.zero.obs~parameters.name,nrow=3, #perc.zero.obs
    scales="free") +
  labs(y="Propability to recover parameter",
       x="",
       fill="functional form") +
  scale_fill_colorblind(labels=c("1.Traditional",
                                 "2.Linear",
                                 "3.Exponent",
                                 "4.Sigmoid"))  +
  theme_bw()

ggsave(ProptoRecoverParam,
       files="ProptoRecoverParam.pdf")
ggsave(ProptoRecoverParam,
       files="ProptoRecoverParam_n.obs.pdf")
ggsave(ProptoRecoverParam,
       files="ProptoRecoverParam_perc.zero.obs.pdf")

#---- Difference between parameter and set value---
diff.df <- df.ite.param%>%
  gather(c("lambda",
           "alpha_init_i", "alpha_init_j",
           "alpha_slope_i","alpha_slope_j",
           "alpha_c_i","alpha_c_j",
           "N_opt_i","N_opt_j"), key="parameter",value="estimate") %>%
  left_join(params_torec.df %>%
              gather(c("lambda",
                       "alpha_init_i", "alpha_init_j",
                       "alpha_slope_i","alpha_slope_j",
                       "alpha_c_i","alpha_c_j",
                       "N_opt_i","N_opt_j"), key="parameter",value="torevover")) %>%
  mutate(diff = estimate- torevover)

Diff.plot <- diff.df %>%
  filter((function.int==1 & parameter %in%param.fctlist[[1]])|
           (function.int==2 & parameter %in% param.fctlist[[2]])|
           (function.int==3 & parameter %in%param.fctlist[[3]])|
           (function.int==4 & parameter %in%param.fctlist[[4]])) %>%
  mutate(parameter = factor(parameter,
                            levels=c("lambda",
                                     "alpha_init_i", "alpha_init_j",
                                     "alpha_slope_i","alpha_slope_j",
                                     "alpha_c_i","alpha_c_j",
                                     "N_opt_i","N_opt_j"))) %>%
  ggplot(aes(group=function.int,x=diff )) +
  geom_density_ridges(aes(y=as.factor(function.int),
                          fill=as.factor(function.int)))+ 
  facet_wrap( #.~parameter,
    #n.obs~parameter,nrow=4, #perc.zero.obs
    perc.zero.obs~parameter,nrow=3, #perc.zero.obs
    scales="free") +
  guides(fill="none") +
  scale_fill_colorblind() + 
  labs(x="Difference between estimated parameter and set parameter",
       y="functional form") +
  theme_bw()

Diff.plot

ggsave(Diff.plot,
       files="Diff.plot.pdf")
ggsave(Diff.plot,
       files="Diff.plot_n.obs.pdf")
ggsave(Diff.plot,
       files="Diff.plot_perc.zero.obs.pdf")
