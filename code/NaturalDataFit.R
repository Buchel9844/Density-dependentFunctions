#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 0. METADATE: Create meta data file to Summarise our data----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 0.1. CREATE TEMPLATES
create_spice() # creates template from directory

# 0.2. EDIT TEMPLATES
#attributes.csv: This is where most of the user data entry will take place. 
# For each variable, its name, units, and a written description are filled in.
edit_attributes()

# access.csv: Includes a row for each file that was read in, 
# and documents the name of each file and its format.
edit_access()

# creators.csv: One row for each creator, and gives their affiliation,
# contact email, ORCID, etc
edit_creators() #

#biblio.csv: Citation information about the project, 
# as much or as little data as possible can be included, 
# but if things like bounding box coordinates are not included, 
# then when the website is generated there will not be a bounding box map generated
edit_biblio()


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
library(car)
library(loo)
library(wesanderson) # for color palette
library(ggthemes) 
#---- 1.2. Import the competitive data ----
competition <- read.csv("/Users/lisabuche/Documents/Projects/Bendering/data/Bendering2022_competition.csv",
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
fecundity.df <- read.csv("/Users/lisabuche/Documents/Projects/Bendering/data/Bendering2022_seeds.csv",
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
# tot number of flowers on individual does not change per capita fecundity - or not enough data to elaborate that

fecundity.summary <- fecundity.df %>%  
  group_by(year, focal) %>%  
  summarise(mean.seed=mean(viable.seeds), st.dev.seed = sd(viable.seeds),
          n= length(viable.seeds),
          mean.flower=mean(tot.flowers),
          st.dev.seed= sd(tot.flowers))%>%  
  ungroup()

focal.levels <- levels(as.factor(competition$focal))
competition.seeds = data.frame()

for(focal in focal.levels){
  competition.n <- subset(competition, competition$focal == focal) 
competition.n$seeds <-  competition.n$flower*
  rnorm(nrow(competition.n),
        mean=fecundity.summary$mean.seed[which(fecundity.summary$focal==focal)],
        sd = fecundity.summary$st.dev.seed[which(fecundity.summary$focal==focal)])
competition.seeds <- rbind(competition.seeds,competition.n)
}
ggplot(competition.seeds,aes(x=seeds, color= focal) ) + geom_density()

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
  #     paste0("results/stan/Competition_seeds_2022_",focal,".csv"))
  Spcompetition[[focal]] <- competition.n
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2. Run the model for each focal----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run.diagnostic =1
run.stan = 1
for(Code.focal in focal.levels){
for (function.int in c(1:4)){
  print(paste(Code.focal,", function",function.int))

    function.vec <- c(0,0,0,0)
  function.vec[function.int] <- 1
 
 # data for the focal
  SpDataFocal <- Spcompetition[[Code.focal]]
 
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
  NamesSpToKeep <-names(SpToKeep)
  S <- sum(SpToKeep)
  #SpMatrix <- matrix(NA, nrow = N, ncol = S)
  #i <- 1
  #for(s in 1:ncol(AllSpAbunds)){
    #if(SpToKeep[s] == 1){
     # SpMatrix[,i] <- AllSpAbunds[,s]
    #  i <- i + 1
   # }else{next}
  #}
  S <- 1+length(focal.levels)
  SpMatrix <- as.data.frame(matrix(NA, nrow = N, ncol = S ))
  #i <- 1
  for(i in 1:N){
   SpMatrix[i,1] <- sum(AllSpAbunds[i,NamesSpToKeep[!NamesSpToKeep %in% focal.levels]])
   SpMatrix[i,c(2:S)] <- AllSpAbunds[i,focal.levels]
  }
  names(SpMatrix) <-c("Neighbours",focal.levels)
  SpMatrix <- as.matrix(SpMatrix)
    #SpMatrix <-round((SpMatrix/max(SpMatrix))*100) #scale all the interaction between 0 and 100
  #if(max(SpMatrix) == 100){print("scale SpMatrix_plant correct")}
  
  SpNames <- AllSpNames[SpToKeep]
  SpNames <-c("Neighbours",focal.levels)
  #assign(paste0("SpNames_",FocalPrefix),
  #     SpNames)
  Intra <- ifelse(SpNames == Code.focal, 1, 0)
  
  # max fecundity
  Nmax <- c(SpMatrix[which.max(Fecundity),])

  # Upper bound intrinsic fecundity
  U <- ceiling(log(fecundity.summary$mean.seed[which(fecundity.summary$focal==Code.focal)]))
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
  
  
  DataVec <- list(N=N, 
                  S=S,
                  U= U,
                  Nmax=Nmax,
                  Fecundity=Fecundity,
                  SpMatrix =SpMatrix,
                  Intra=Intra,
                  run_estimation=run_estimation,
                  alphaFunct1=alphaFunct1,
                  alphaFunct2=alphaFunct2,
                  alphaFunct3=alphaFunct3,
                  alphaFunct4=alphaFunct4
                  )
  
  ##---- 3.2. Run  final fit ----
  # Now run a fianl fit of the model to assess parameter 
  print("Natural Fit beginning")
  
  #install.packages("codetools")
  library("codetools")
  options(mc.cores = parallel::detectCores())
  # defiining initial value of lambda to help with computation failure
  # divide by the maximum expected for lambda = U 
  list.init <- function(...)list(lambdas = array(as.numeric(rnorm(1,
                                                                  mean=log(fecundity.summary$mean.seed[which(fecundity.summary$focal==Code.focal)]),
                                                                  sd = abs(log(fecundity.summary$st.dev.seed[which(fecundity.summary$focal==Code.focal)])))/
                                                              U), 
                                                            dim = 1))

  if( run.stan == 1){                               
  FinalFit <- stan(file = "code/DensityFunct_Final.stan", 
                   data = DataVec,
                   init=  list.init,
                   warmup= 500,
                   iter = 1000, 
                   init_r = 2,
                   chains = 3,
                   control=list(max_treedepth=15),
                   seed= 165)
  
  
  save(file= paste0("results/stan/FinalFit_",Code.focal,"_",alpha.function,".rds"),
       FinalFit)
  }
  #load(paste0("results/stan/FinalFit_",Code.focal,"_",alpha.function,".rds"))
  
 load(paste0("results/stan/FinalFit_",Code.focal,"_",alpha.function,".rds"))

 FinalPosteriors <- rstan::extract(FinalFit)
  
  print("Final Fit done")
  
  #---- 3.3. Final fit posterior check and behavior checks---- 
  if( run.diagnostic == 1){  
  ##### Diagnostic plots and post prediction 
  pdf(paste0("figures/stan/FinalFit_",Code.focal,"_",alpha.function,".pdf"))
  # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
  source("code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
  # check the distribution of Rhats and effective sample sizes 
  ##### Posterior check
  stan_post_pred_check(FinalPosteriors,"F_hat",Fecundity,
                       paste0("results/stan/PostFec_",Code.focal,"_",alpha.function,".csv.gz")) 
  
  log_lik_2 <- loo::extract_log_lik(FinalFit, 
                                    parameter_name = "F_sim", 
                                    merge_chains = F)
  
  r_eff <- loo::relative_eff(exp(log_lik_2), cores = 2) 
  
  loo_1 <- loo::loo(log_lik_2 , r_eff = r_eff, cores = 2)
  
  print(loo_1)
  
  # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
  hist(summary(FinalFit)$summary[,"Rhat"],
       main = paste("Finat Fit: Histogram of Rhat for",
                    Code.focal," and ",alpha.function))
  hist(summary(FinalFit)$summary[,"n_eff"],
       main = paste("Finat Fit: Histogram of Neff for",
                    Code.focal," and ",alpha.function))
  
  # plot the corresponding graphs
  trace <- stan_trace(FinalFit, pars=c('lambdas','c','alpha_initial','alpha_slope','c'),
                      inc_warmup = TRUE)
  print(trace)
  dens <- stan_dens(FinalFit, 
                    pars=c('lambdas','c','alpha_initial','alpha_slope','c'))
  print(dens)
  splot <- stan_plot(FinalFit, 
                     pars=c('lambdas','c','alpha_initial','alpha_slope','c'))
  print(splot)
  sampler_params <- get_sampler_params(FinalFit, inc_warmup = TRUE)
  summary(do.call(rbind, sampler_params), digits = 2)
  pairs(FinalFit, pars = c("lambdas",'alpha_initial',
                           'alpha_slope','c'))
  
  # Next check the correlation among key model parameters and identify any
  #pairs(FinalFit, pars = c("lambdas",'alpha_initial','alpha_slope','c'))
  
  
  dev.off()
  }
  #---- 3.3. Extract coefficients ----
  
  assign(paste0("Parameters_",Code.focal,"_",alpha.function),
       list(DataVec = DataVec,
       alpha_value=FinalPosteriors$alpha_value,      
       lambda =  FinalPosteriors$lambda_ei,
       alpha_slope= FinalPosteriors$alpha_slope,
       alpha_init =FinalPosteriors$alpha_init,
       alpha_c = FinalPosteriors$c
       ))
  
  save(list =paste0("Parameters_",Code.focal,"_",alpha.function),
       file = paste0("results/stan/Parameters_",Code.focal,"_",alpha.function,".RData"))
  
  #---- 3.3. Extraction fecundity---
  
  assign(paste0("Fsim_",Code.focal,"_",alpha.function),
         list(DataVec = DataVec,
              F_sim=FinalPosteriors$F_sim
         ))
  
  save(list =paste0("Fsim_",Code.focal,"_",alpha.function),
       file = paste0("results/stan/Fsim_",Code.focal,"_",alpha.function,".RData"))

}
}

#---- 4. Figures ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("results/stan/NaturalData_PostFecunditydistribution.csv.gz")
Post_Fecunditydistribution <- PostFecunditydistribution %>%
  filter(!is.na(Fec))

#---- 4.0. check models ----
ModelCheck <- NULL
for(Code.focal in focal.levels){
    for (function.int in c(1:4)){ # c(1:4)
      if (Code.focal == "TROR") next
      print(paste(Code.focal,", function",function.int))
      
      alpha.function <- paste0("function_",function.int)
      
      load(paste0("results/stan/FinalFit_",Code.focal,"_",alpha.function,".rds"))
      
      mc <- data.frame(focal =Code.focal, 
                        function.int = function.int,
                        Rhat = max(summary(FinalFit)$summary[,"Rhat"],na.rm =T),
                        Neff = min(summary(FinalFit)$summary[,"n_eff"],na.rm = T))
      ModelCheck <- bind_rows(ModelCheck,mc)
      remove(FinalFit)
    }
  }
View(ModelCheck)


#---- 4.1 Alpha distribution ----

neighbours.vec <- c("neighbours",focal.levels)

df_alpha_all <- data.frame(focal=NA)
for(Code.focal in focal.levels){
  for (function.int in c(1:4)){
    if (Code.focal == "TROR") next
    df_alpha <- data.frame()
    function.vec <- c(0,0,0,0)
    function.vec[function.int] <- 1
    alpha.function <- paste0("function_",which(function.vec==1))

    load(paste0("results/stan/Parameters_",Code.focal,"_",alpha.function,".RData"))
    parameters <- get(paste0("Parameters_",Code.focal,"_",alpha.function))
    N <- parameters[["DataVec"]]$N
    Nmax <- parameters[["DataVec"]]$Nmax
    neighbours.vec <- c("neighbours",focal.levels)
    df_alpha_init <- data.frame(observation = c(1:N),
                                 focal = Code.focal,
                                 function.int = function.int)
     df_alpha_init <- NULL
     for(n in 1:length(neighbours.vec)){
       df_alpha_init_n <- data.frame(focal = Code.focal,
                                     Nmax=Nmax[n],
                                     function.int = function.int,
                                     number = c(1:length(parameters[["alpha_init"]][,n])),
                                     neigh = neighbours.vec[n],
                                     alpha_init= parameters[["alpha_init"]][,n])
       df_alpha_init <- bind_rows(df_alpha_init,df_alpha_init_n)
        }

    if (function.int > 1){
      df_alpha_slope <- NULL
      for(n in 1:length(neighbours.vec)){
        df_alpha_slope_n <- data.frame(focal = Code.focal,
                                      neigh = neighbours.vec[n],
                                      function.int = function.int,
                                      number = c(1:length(parameters[["alpha_slope"]][,n])),
                                      alpha_slope = parameters[["alpha_slope"]][,n])
        df_alpha_slope <- bind_rows(df_alpha_slope,df_alpha_slope_n)
      }
      
     if(function.int > 2){
       df_alpha_c <- NULL
       for(n in 1:length(neighbours.vec)){
         df_alpha_c_n <- data.frame(focal = Code.focal,
                                        neigh = neighbours.vec[n],
                                        function.int = function.int,
                                        number = c(1:length(parameters[["alpha_c"]][,n])),
                                        alpha_c = parameters[["alpha_c"]][,n])
         df_alpha_c <- bind_rows(df_alpha_c,df_alpha_c_n)
       }
    
    df_alpha <- full_join(df_alpha_c,full_join(df_alpha_init,df_alpha_slope))
     }else{
    df_alpha <- full_join(df_alpha_init,df_alpha_slope)
     }
    }else{
    df_alpha <- df_alpha_init
    }
     df_alpha_all <- full_join(df_alpha_all,df_alpha)
  }
}
df_alpha_all <- df_alpha_all[-1,] # nrow = 120 == 4(fct) * 5(species) * 6(neighbours)
str(df_alpha_all)




df_funct_alpha <- NULL
for (n in 1:nrow(df_alpha_all)){
  print(n)
  df_funct_alpha_n <- df_alpha_all[n,]
  df_funct_alpha_n <- data.frame(density=c(0:10),df_funct_alpha_n)
  if(df_funct_alpha_n$function.int[1]==1){
    df_funct_alpha_n$alpha_value <- df_funct_alpha_n$alpha_init
    
    }
  if(df_funct_alpha_n$function.int[1]==2){
    df_funct_alpha_n$alpha_value <- alpha_function2(df_funct_alpha_n$alpha_init,
                                                  df_funct_alpha_n$alpha_slope,
                                                  df_funct_alpha_n$density,
                                                  df_funct_alpha_n$Nmax)
    
  }
  if(df_funct_alpha_n$function.int[1]==3){
    df_funct_alpha_n$alpha_value <- alpha_function3(df_funct_alpha_n$alpha_init,
                                                  df_funct_alpha_n$alpha_slope,
                                                  df_funct_alpha_n$alpha_c,
                                                  df_funct_alpha_n$density,
                                                  df_funct_alpha_n$Nmax)
    
  }
  if(df_funct_alpha_n$function.int[1]==4){
    df_funct_alpha_n$alpha_value <- alpha_function4(df_funct_alpha_n$alpha_init,
                                                  df_funct_alpha_n$alpha_slope,
                                                  df_funct_alpha_n$alpha_c,
                                                  df_funct_alpha_n$density,
                                                  df_funct_alpha_n$Nmax)
  }
  
  df_funct_alpha <- bind_rows(df_funct_alpha,df_funct_alpha_n )
  
}

write.csv(df_funct_alpha,
          "results/df_funct_alpha.csv.gz")


cbp2 <- c("#000000", "#E69F00", "#56B4E9","#F0E442","#009E73",
          "#CC79A7", "#0072B2", "#D55E00")

cb2_focal <- c("#E69F00", "#56B4E9","#009E73","#CC79A7")

plot.list.alpha <- list()



# percentage of positive interaction with function 4
100*nrow(df_funct_alpha[which(df_funct_alpha$function.int == 4 &
                       df_funct_alpha$alpha_value > 0),])/
  nrow(df_funct_alpha[which(df_funct_alpha$function.int == 4),])


for( n in 1:4){
  df <- df_funct_alpha[which(df_funct_alpha$focal == focal.levels[n]),]
  plot.list.alpha[[n]] <- ggplot(df,aes(x=density, y= alpha_value,
                                      color=neigh,fill=neigh)) +
    stat_smooth(method = 'gam',se = TRUE,level =0.95) +
    facet_wrap(. ~ function.int, strip.position = "top",
               nrow=1, ncol=4, scales = "fixed") +
    scale_color_manual("neighbours identity",values=cbp2[2:7]) + 
    scale_fill_manual("neighbours identity",values=cbp2[2:7]) + 
    theme_bw() +
    labs(title = element_text(focal.levels[n],
                              color=cb2_focal[n]),
         y= "Effect of neighbours on focal",
         x="density of neighbours") +
    rremove("ylab") + rremove("xlab") + 
    geom_hline(yintercept=0, color="black", linetype="dashed") +
    theme(plot.title = element_text(color=cb2_focal[n],
                                       vjust = -5),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text =element_text(color="black"),
          legend.key= element_rect(fill = "white"),
          plot.margin = unit(c(-1,0.2,0,0.2), 'lines')) 
}


plot.alpha_all <- ggarrange(plotlist = plot.list.alpha,
                            labels = NULL,
                            align = "hv",
                            font.label = list(size = 10, color = "black", 
                                              face = "bold", family = NULL, position = "top"),
                            nrow=4, common.legend=T, legend="right")
library(grid)
plot.alpha_all <- annotate_figure(plot.alpha_all, left = textGrob("Effect of neighbours on focal", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("density of neighbours", gp = gpar(cex = 1.3)))
ggsave(plot.alpha_all,
       file = "figures/NatData_plot_alpha_all.pdf")

plot.list.alpha <- list()

for( n in 1:4){
  df <- df_funct_alpha[which(df_funct_alpha$focal == focal.levels[n] &
                               (df_funct_alpha$function.int == 1|
                                  df_funct_alpha$function.int == 4)),]
  plot.list.alpha[[n]] <- ggplot(df,aes(x=density, y= alpha_value,
                                        color=neigh,fill=neigh)) +
    stat_smooth(method = 'gam',se = TRUE,level =0.95) +
    facet_wrap(. ~ function.int, strip.position = "top",
               nrow=1, ncol=2, scales = "fixed") +
    scale_color_manual("neighbours identity",values=cbp2[2:7]) + 
    scale_fill_manual("neighbours identity",values=cbp2[2:7]) + 
    theme_bw() +
    labs(title = element_text(focal.levels[n],
                              color=cb2_focal[n]),
         y= "Effect of neighbours on focal",
         x="density of neighbours") +
    rremove("ylab") + rremove("xlab") + 
    geom_hline(yintercept=0, color="black", linetype="dashed") +
    theme(plot.title = element_text(color=cb2_focal[n],
                                    vjust = -5),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text =element_text(color="black"),
          legend.key= element_rect(fill = "white"),
          plot.margin = unit(c(-1,0.2,0,0.2), 'lines')) 
}



plot.alpha_short <- ggarrange(plotlist = plot.list.alpha,
                            labels = NULL,
                            align = "hv",
                            font.label = list(size = 10, color = "black", 
                                              face = "bold", family = NULL, position = "top"),
                            nrow=4, common.legend=T, legend="right")
library(grid)
plot.alpha_short <- annotate_figure(plot.alpha_short, left = textGrob("Effect of neighbours on focal", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                  bottom = textGrob("density of neighbours", gp = gpar(cex = 1.3)))
ggsave(plot.alpha_short,
       file = "figures/NatData_plot_alpha_short.pdf")


#---- 4.2. Population Projections ----

df_param_all <- data.frame(focal=NA)
for(Code.focal in focal.levels){
  for (function.int in c(1:4)){
    if (Code.focal == "TROR") next
    df_param <- data.frame()
    function.vec <- c(0,0,0,0)
    function.vec[function.int] <- 1
    alpha.function <- paste0("function_",which(function.vec==1))
    
    load(paste0("results/stan/Parameters_",Code.focal,"_",alpha.function,".RData"))
    parameters <- get(paste0("Parameters_",Code.focal,"_",alpha.function))
    N <- parameters[["DataVec"]]$N
    Nmax <- parameters[["DataVec"]]$Nmax
    neighbours.vec <- c("neighbours",focal.levels)
    df_param_init <- data.frame(observation = c(1:N),
                                focal = Code.focal,
                                function.int = function.int)
    df_param_init <- NULL
    for(n in 1:length(neighbours.vec)){
      df_param_init_n <- data.frame(focal = Code.focal,
                                    Nmax=Nmax[n],
                                    function.int = function.int,
                                    neigh = neighbours.vec[n],
                                    lambda= mean(parameters[["lambda"]]),
                                    alpha_init= mean(parameters[["alpha_init"]][,n]))
      df_param_init <- bind_rows(df_param_init,df_param_init_n)
    }
    
    if (function.int > 1){
      df_param_slope <- NULL
      for(n in 1:length(neighbours.vec)){
        df_param_slope_n <- data.frame(focal = Code.focal,
                                       neigh = neighbours.vec[n],
                                       function.int = function.int,
                                       alpha_slope = mean(parameters[["alpha_slope"]][,n]))
        df_param_slope <- bind_rows(df_param_slope,df_param_slope_n)
      }
      
      if(function.int > 2){
        df_param_c <- NULL
        for(n in 1:length(neighbours.vec)){
          df_param_c_n <- data.frame(focal = Code.focal,
                                     neigh = neighbours.vec[n],
                                     function.int = function.int,
                                     alpha_c = mean(parameters[["alpha_c"]][,n]))
          df_param_c <- bind_rows(df_param_c,df_param_c_n)
        }
        
        df_param <- full_join(df_param_c,full_join(df_param_init,df_param_slope))
      }else{
        df_param <- full_join(df_param_init,df_param_slope)
      }
    }else{
      df_param <- df_param_init
    }
    df_param_all <- full_join(df_param_all,df_param)
  }
}
df_param_all <- df_param_all[-1,] # nrow = 120 == 4(fct) * 5(species) * 6(neighbours)
str(df_param_all)
df_projection <- NULL

levels(as.factor(df_param_all$lambda))

for (function.int in c(1:4)){
  param.df <- df_param_all[which(df_param_all$function.int==function.int),]
  param.df$g = 1
  param.df$s = 1
  state = c(1)
  
  df_projection_n <- Ricker_solution_NatData(gens = 250,
                                             state,pars = param.df)
  df_projection_n$function.int <- function.int
  df_projection <-  bind_rows(df_projection, df_projection_n)
}

plot_projection <- df_projection %>%
  gather(all_of(c(focal.levels,"neighbours")), key="species", value="abundance") %>%
  ggplot(aes(y=abundance, x= time, group=species, color=species)) +
  geom_smooth() + scale_y_log10() +
  facet_grid(.~as.factor(function.int),
             scales="free") +
  theme_bw()

ggsave(plot_projection, 
       file = "figures/NatData_plot_projection.pdf")
#---- 4.3. Compare model for each focal  ----
# reference : Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing. 27(5), 1413–1432. :10.1007/s11222-016-9696-4. online, arXiv preprint arXiv:1507.04544.
model.loo <- list()

for(Code.focal in focal.levels){
  for (function.int in 1:4){ # c(1:4)
    if (Code.focal == "TROR") next
    
    if(ModelCheck[which(ModelCheck$function.int==function.int &
                       ModelCheck$focal==Code.focal),"Rhat"] < 0 | 
      ModelCheck[which(ModelCheck$function.int==function.int &
                         ModelCheck$focal==Code.focal),"Neff"] < 100) next
      
      print(paste(Code.focal,", function",function.int))
      
      load(paste0("results/stan/FinalFit_",Code.focal,"_function_",function.int,".rds"))
      # Extract pointwise log-likelihood
      # using merge_chains=FALSE returns an array, which is easier to 
      # use with relative_eff()
      log_lik <- loo::extract_log_lik(FinalFit, 
                                      parameter_name = "F_sim", 
                                      merge_chains = F)
      #as of loo v2.0.0 we can optionally provide relative effective sample sizes
      # when calling loo, which allows for better estimates of the PSIS effective
      # sample sizes and Monte Carlo error
      r_eff <- loo::relative_eff(log_lik, cores = 2) 
      # preferably use more than 2 cores (as many cores as possible)
      # will use value of 'mc.cores' option if cores is not specified
      model.loo[[paste0(Code.focal,"_function_",function.int)]] <- loo::loo(log_lik,
                                                                            r_eff = r_eff, 
                                                                            cores = 2)
      remove(FinalFit)
  }
}
Se_loo_model <- NULL
for(Code.focal in focal.levels){
  if (Code.focal == "TROR") next
    comp <- loo:loo_compare(model.loo[[paste0(Code.focal,"_function_1")]], 
                        model.loo[[paste0(Code.focal,"_","function_2")]],
                        model.loo[[paste0(Code.focal,"_","function_3")]],
                        model.loo[[paste0(Code.focal,"_","function_4")]])
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
  save(Se_loo_model,
       file = paste0("results/model.loo.",scenario,".RData"))
  
write.csv(Se_loo_model,
          paste0("results/Se_loo_model.csv"))

#---- 4.3. Figure of fecundity estimation distribution ----

  pdf(paste0("figures/NatData_PostFecundity_distribution.pdf")) 
  par(oma=c(0,0,0,0), mar=c(2.25,2.5,2,1))
  layout(mat = matrix(
    1:16,
    byrow=TRUE,
    nrow = 4,
    ncol = 4
  ),
  heights = rep(4,8),
  widths = rep(3,8)
  )
  
  for(Code.focal in focal.levels){ #,"j"
    for (function.int in c(1:4)){ # c(1:4)
      if (Code.focal == "TROR") next
     #col2 <- c("black","#CC79A7","#E69F00","#009E73") # function 
      col2 <- wes_palette("FantasticFox1", n = 5)
      #col1 <- c("blue","red") # focal species
      # loo values
      #value.se <- round(Se_loo_model$se_diff[which(Se_loo_model$model == function.int &
      #                                               Se_loo_model$focal == Code.focal)],digits=2)
      # fecundity obs
      load(paste0("results/stan/Fsim_",Code.focal,"_function_",function.int,".RData"))
      # fecundity generated
      Fec_df <- data.frame(Obs = 1:get(paste0("Fsim_",Code.focal,"_function_",function.int))[["DataVec"]]$N,
                           Fec =get(paste0("Fsim_",Code.focal,"_function_",function.int))[["DataVec"]]$Fecundity)
      
      
      load(paste0("results/stan/FinalFit_",Code.focal,"_function_",function.int,".rds"))
      FinalPosteriors <- rstan::extract(FinalFit)
      stan_post_pred_check_all(FinalPosteriors,"F_hat",
                               Fec_df$Fec,
                               paste0("Species ",Code.focal,",function ",function.int),
                               "black",
                               col2[1+function.int]#,#value.se
                               ) 
      
      
    }
  }
  dev.off()



