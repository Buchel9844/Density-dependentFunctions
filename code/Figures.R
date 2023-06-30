# This script will make the figures of chap 2, number- XXX of the manuscript

library(ggplot2)
library(ggridges)
library(tidyverse)
##########################################################################################################
# 1. Simulation data
##########################################################################################################
#---- 0. check models ----
ModelCheck <- NULL
for( scenario in c("low","medium","high")){
  for(Code.focal in c("i","j")){ #,"j"
    for (function.int in c(1:4)){ # c(1:4)
      
      print(paste(scenario, Code.focal,", function",function.int))
    
      load(paste0("results/FinalFit_",scenario,"_",Code.focal,"_function_",function.int,".rds"))
      mc <- data.frame( scenario = scenario, focal =Code.focal, 
                        function.int = function.int,
                        Rhat = max(summary(FinalFit)$summary[,"Rhat"],na.rm =T),
                        Neff = min(summary(FinalFit)$summary[,"n_eff"],na.rm = T))
      ModelCheck <- bind_rows(ModelCheck,mc)
      remove(FinalFit)
    }
  }
}
  View(ModelCheck)
#---- 1.0. Compare simulation with estimate alpha for ricker model  ----
load("results/Alphadistribution.neighbours.Ricker.csv.gz")

sim_alpha_specific <- as.data.frame(simul_data_ricker$sim_alpha_specific) %>%
  rownames_to_column(var = "focal") %>%
  gather(alphai,alphaj,key="neighbours" ,value="alpha_sim") %>%
  mutate(focal = case_when(focal == "focali" ~ "species i",
                           focal == "focalj" ~ "species j"),
         neighbours = case_when(neighbours == "alphai" ~ "species i",
                                neighbours == "alphaj" ~ "species j"))

write.csv(full_join(sim_alpha_specific ,
                    Alphadistribution.neighbours[which(Alphadistribution.neighbours$density.function=="function_2" & 
                                                         Alphadistribution.neighbours$abundance.neighbours == 1),]),
          file= "results/Alphadistribution_RickerModel.csv"
)

AlphadistributionGraph.Ricker <- ggplot() + 
  geom_smooth(data=Alphadistribution.neighbours.Ricker, aes(x= abundance.neighbours,
                                                       y=alpha_mean,
    color=density.function,
    fill=density.function),
    alpha=0.5,se=F) +
  geom_point(data=Alphadistribution.neighbours.Ricker,aes(
    x= abundance.neighbours,
    y=alpha_mean,
    color=density.function,
    fill=density.function)) +
  #ylim(-0.1,0.1)+
  geom_hline(yintercept=0,color="dark grey") +
  geom_hline( data=sim_alpha_specific ,aes(yintercept=alpha_sim),color="blue") + 
  geom_errorbar(data=Alphadistribution.neighbours.Ricker, 
                aes(x= abundance.neighbours,
                    y=alpha_mean,
                    ymin=alpha_mean-alpha_sd^2, ymax=alpha_mean+alpha_sd^2,
                    color=density.function), width=.2,
                position=position_dodge(0.05)) +
  facet_grid( focal ~ neighbours,scale="free", switch="both") + theme_bw() + 
  guides(fill="none") +
  scale_color_manual("Density-dependent functions",values=c("black","#CC79A7","#E69F00","#009E73")) +
  labs(title="Direct interactions distributions of the 4 density-dependent functions from RICKER simulation data", 
       y="Resulting effect", x=" Neighbour density ")

ggsave("figures/AlphadistributionGraph.Ricker.pdf",
       plot = AlphadistributionGraph.Ricker
)


#---- 1.1. Figure of density-dependent function ----
for( scenario in c("low","medium","high")){
     
load(paste0("results/Alphadistribution.",scenario,".neighbours.csv.gz"))

Alphadistribution.neighbours$interaction <- paste0("alpha_",
                                                     substr(Alphadistribution.neighbours$focal,9,9),
                                                     substr(Alphadistribution.neighbours$neighbours,9,9))

AlphadistributionGraph <- ggplot(Alphadistribution.neighbours, aes(x= c(abundance.neighbours - Nmax),
                                                           y=alpha_mean)) + 
  geom_hline(yintercept=0,color="dark grey") +
  stat_smooth(geom='line', 
              aes(
                color=density.function),
              alpha=0.7,se=F,size=2)+
  geom_point(aes(
    color=density.function,
    fill=density.function),alpha=0.8,size=2) +
  #ylim(-0.1,0.1)+
  #xlim(0,max(Alphadistribution.neighbours$abundance.neighbours)) + 
  geom_errorbar(aes(ymin=alpha_mean-alpha_sd^2, ymax=alpha_mean+alpha_sd^2,
                    color=density.function), width=.5) +
  guides(fill="none") +
  scale_color_manual("Density-dependent functions",values=c("black","#CC79A7","#E69F00","#009E73")) +
  labs(title="Direct interactions distributions of the 4 density-dependent functions", 
       y="Resulting effect", x=" Neighbour density ") + 
  geom_text(mapping = aes(x = max(abundance.neighbours - Nmax), 
                          y = max(alpha_mean)-0.1, label = interaction),
            hjust   = 1) +
  facet_grid( focal ~ neighbours,scale="free", switch="both") + theme_bw() 
  

ggsave(paste0("figures/AlphadistributionGraph",scenario,".pdf"),
       plot = AlphadistributionGraph
)

}


#---- 1.2. Figure of fecundity estimation distribution ----
for( scenario in c("low","medium","high")){
  
load( paste0("results/Fecunditydistribution.",scenario,".csv.gz"))

load( paste0("results/PostFecunditydistribution.",scenario,".csv.gz"))

simdata <- read.csv(paste0("results/Generate.simulated.data.",scenario,".csv"))
simdata <- simdata[which(simdata$time <= time.exp),]

Fecunditydistribution.n <- data.frame(seed=simdata$fecundity,
                                    focal=simdata$focal,
                                    obervation= NA )

Post_Fecunditydistribution <- PostFecunditydistribution %>%
  filter(!is.na(Fec))

PostFecundityGraph <- ggplot(PostFecunditydistribution) +
  geom_density(aes(x=Fec,
                  group= obs,
                  color=as.factor(alpha.function)),
               alpha=0.2) +
  geom_density(data = Fecunditydistribution.n , aes(x=seed, 
                                          color=as.factor(focal)),
               alpha=0.5) +
  facet_wrap(focal ~  alpha.function, nrow=2, scales = "free" ) +
  scale_x_continuous(limits = function(x) {x + c(-10, 10)}) +
  labs(title = paste0(" Fecundity distributions for predictions compared to initial data for ",scenario), y ="Density",
       x="Fecundity") +  theme_bw() + guides(color="none") + 
  scale_color_manual(values=c("black","#CC79A7","#E69F00","#009E73",
                              "blue","red"))

ggsave(paste0("figures/PostFecundityGraph",scenario,".pdf"),
       plot = PostFecundityGraph)
remove(PostFecunditydistribution)
remove(simdata)
}
#---- 1.2.bis. Figure of fecundity estimation distribution ----
for( scenario in c("low","medium","high")){
pdf(paste0("figures/PostFecundity_distribution_",scenario,".pdf")) 
  par(oma=c(0,0,0,0), mar=c(2.25,2.5,2,1))
  layout(mat = matrix(
    1:8,
    byrow=TRUE,
    nrow = 2,
    ncol = 4
  ),
  heights = c(6,6),
  widths = rep(3,8)
  )
  
  for(Code.focal in c("i","j")){ #,"j"
    for (function.int in c(1:4)){ # c(1:4)
      
load(paste0("results/FinalFit_",scenario,"_",Code.focal,"_function_",function.int,".rds"))
simdata <- read.csv(paste0("results/Generate.simulated.data.",scenario,".csv"))
simdata <- simdata[which(simdata$time >= time.exp.min & 
                           simdata$time < time.exp.max &
                           simdata$focal == Code.focal),]

col2 <- c("black","#CC79A7","#E69F00","#009E73")
col1 <- c("blue","red")
if(Code.focal=="i"){int.focal <- 1}else{int.focal <- 2}    
value.se <- round(Se_loo_model$se_diff[which(Se_loo_model$model == function.int &
                                         Se_loo_model$scenario == scenario &
                                         Se_loo_model$focal == Code.focal)],digits=2)

FinalPosteriors <- rstan::extract(FinalFit)
stan_post_pred_check_all(FinalPosteriors,"F_hat",
                     simdata$fecundity[!is.na(simdata$fecundity)],
                     paste0("Species ",Code.focal,",function ",function.int),
                     col1[int.focal],
                     col2[function.int ],
                     value.se) 


    }
  }
  dev.off()
}
  
  
#---- 2.0. Compare model for each focal  ----
# reference : Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing. 27(5), 1413–1432. :10.1007/s11222-016-9696-4. online, arXiv preprint arXiv:1507.04544.
model.loo <- list()
for( scenario in c("low","medium","high")){
for(Code.focal in c("i","j")){ #,"j"
  for (function.int in c(1:4)){ # c(1:4)
    
load(paste0("results/FinalFit_",scenario,"_",Code.focal,"_function_",function.int,".rds"))
    # Extract pointwise log-likelihood
    # using merge_chains=FALSE returns an array, which is easier to 
    # use with relative_eff()
log_lik <- loo::extract_log_lik(FinalFit, 
                                  parameter_name = "F_sim", 
                                  merge_chains = F)
#as of loo v2.0.0 we can optionally provide relative effective sample sizes
# when calling loo, which allows for better estimates of the PSIS effective
# sample sizes and Monte Carlo error
r_eff <- loo::relative_eff(exp(log_lik), cores = 2) 
# preferably use more than 2 cores (as many cores as possible)
# will use value of 'mc.cores' option if cores is not specified
model.loo[[paste0(scenario,"_",Code.focal,"_function_",function.int)]] <- loo::loo(log_lik, 
                                                               r_eff = r_eff, cores = 2)
remove(FinalFit)
    }
  }
}
Se_loo_model <- NULL
for( scenario in c("low","medium","high")){
for(Code.focal in c("i","j")){ #,"j"
comp <- loo_compare(model.loo[[paste0(scenario,"_",Code.focal,"_","function_1")]], 
                    model.loo[[paste0(scenario,"_",Code.focal,"_","function_2")]],
                    model.loo[[paste0(scenario,"_",Code.focal,"_","function_3")]],
                    model.loo[[paste0(scenario,"_",Code.focal,"_","function_4")]])
model.loo[[Code.focal]] <- comp # The first column shows the difference in ELPD relative to the model with the largest ELPD.

  df_loo <- data.frame(se_diff = model.loo[[Code.focal]][,"se_diff"]) %>%
 rownames_to_column(var = "model") %>%
    mutate(model = stringi::stri_sub(model,from =6,to=6),
           scenario = scenario,
           focal = Code.focal)

  Se_loo_model <- bind_rows(Se_loo_model, df_loo)
}
assign(paste0("model.loo.",scenario),
       model.loo)
list.name <- paste0("model.loo.",scenario)
save(list.name,
     file = paste0("results/model.loo.",scenario,".RData"))
}
write.csv(Se_loo_model,
          paste0("results/Se_loo_model.csv"))


##########################################################################################################
# 2. Empirical data
##########################################################################################################

#---- 2.1. Figure of empirical seed distribution ----

