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
library(grid)
#---- 1.2. Import the competitive data ----

competition <- read.csv("/Users/lisabuche/Documents/Projects/Perenjori/data/Sevenello2022_competition.csv",
                        header = T,stringsAsFactors = F, sep=",",
                        na.strings=c("","NA"))
species.neigh <- names(competition)
species.neigh  <- species.neigh[!species.neigh %in% c("reserve","focal","crop", "site", "plot","treatment","year")]

competition.to.keep <- competition %>%  
  dplyr::select(all_of(c(species.neigh)))%>%
  mutate_at( species.neigh, as.numeric) %>%
  colSums(na.rm=T)
length(competition.to.keep[competition.to.keep == 0]) # check if any is 0


competition <- competition %>% 
  dplyr::filter(treatment =="OP") 

# change the format of the rows to numeric 
competition[species.neigh] <- sapply(competition[species.neigh],as.numeric)

# change na values to 0
competition[is.na(competition)] <- 0

focal.levels <- levels(as.factor(competition$focal))

#---- 1.2. Make Teasorus ----
perenjory_tesaorus <- read.csv("/Users/lisabuche/Documents/Projects/Perenjori/data/perenjory_tesaorus.csv",
                        header = T,stringsAsFactors = F, sep=",",
                        na.strings=c("","NA"))

competition.long.neigh <- competition %>%
  gather(any_of(perenjory_tesaorus$code), key="code",value="abundance") %>%
  left_join(perenjory_tesaorus) 

competition.long.focal <- competition.long.neigh %>%
  filter(focal ==code) %>%
  rename("conspecific"="abundance") %>%
  dplyr::select(-code) 

competition.long <- competition.long.neigh %>%
  aggregate(abundance ~  family + reserve + year + site + focal+ crop + site + plot + treatment, sum )%>%
  spread(family, abundance)  %>%# change this ti have grass /forb or native/exotic
  right_join(competition.long.focal)

# Only run if not family
competition.long <- competition.long.neigh %>%
  aggregate(abundance ~  functional.group + reserve + year + site + focal+ crop + site + plot + treatment, sum )%>%
  spread(functional.group, abundance)  %>%# change this ti have grass /forb or native/exotic
  right_join(competition.long.focal)

#---- 1.4. Import the seeds data ----
fecundity.df <- read.csv("/Users/lisabuche/Documents/Projects/Perenjori/data/Sevenello2022_seeds.csv",
                         header = T,stringsAsFactors = F, sep=",",
                         na.strings=c("","NA"))
head(fecundity.df)
fecundity.df <- fecundity.df %>% 
  dplyr::filter(treatment =="OP")

fecundity.df$reserve <- as.factor(fecundity.df$reserve)

seeds.model <- glm(round(seeds) ~ reserve + focal  ,fecundity.df, 
                     family="poisson")

summary(seeds.model) # random Side effect not signidicant - not considered further
  
seeds.graph <- ggplot(fecundity.df,aes(x=seeds,color=focal)) + 
    geom_density() +
  xlim(0,200)
    
# join both compeetiton and seed 
competition.seeds <- left_join(competition.long,fecundity.df,
                               by=c("site", "focal", "crop", "reserve", "plot", "treatment"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2. Run the model for each focal----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run.diagnostic =1
run.stan = 1
function.grouping = 1
grouping ="family"
for(Code.focal in "LARO"){ #focal.levels
 for (function.int in c(1:4)){
  print(paste(Code.focal,", function",function.int))
    function.vec <- c(0,0,0,0)
  function.vec[function.int] <- 1
  
 # data for the focal
  SpDataFocal <- competition.seeds[which(competition.seeds$focal==Code.focal),] 
  names(SpDataFocal) <- tolower(names(SpDataFocal))
  SpDataFocal <-  SpDataFocal[!is.na(SpDataFocal$seeds),]
  SpDataFocal[is.na(SpDataFocal)] <- 0
  # Next continue to extract the data needed to run the model. 
  N <- as.integer(nrow(SpDataFocal))
  Fecundity <- as.integer(SpDataFocal$seeds)  
  if(grouping =="family"){
  specific.focal <- tolower(SpDataFocal$family[1])
  }else{
    specific.focal <- tolower(SpDataFocal$functional.group[1])
    
  }
  #---- 2. ABUDANCE MATRIX----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #---- 2.1. Interaction (direct) matrix of plant with COMP ----
  
  # Now calculate the total number of plant species to use for the model, discounting
  #       any species columns with 0 abundance. Save a vector of the species names
  #       corresponding to each column for easy matching later.
  AllSpNames <- names(SpDataFocal)[!names(SpDataFocal) %in% c("focal","year","day","month","transect","plot","reserve",
                                                              "block","unique.plot","site","treatment",
                                                              "crop","flower","seeds","individual", "species","origine","functional.group",      
                                                              "family","flowercolor","final.code","flowers",                
                                                              "Seeds")]
  AllSpAbunds <- SpDataFocal %>% 
    dplyr::select(all_of(c(AllSpNames)))%>%
    mutate_at( AllSpNames, as.numeric)
 
   SpTotals <- colSums(AllSpAbunds)
   SpToKeep <- SpTotals > 30
   sum(SpTotals)
   NamesSpToKeep <- names(SpTotals[SpToKeep])
   Stotal <- sum(SpToKeep)
   S <- Stotal
   SpMatrix <- NULL
   SpMatrix <- as_tibble(matrix(NA, nrow = N, ncol = S+1))
   i <- 1
   SpMatrix[,S+1] <- 0
  for(s in 1:length(AllSpNames)){
    if(SpToKeep[s] == 1){
      if(SpTotals[s] < 5){
        SpMatrix[,S+1] <- SpMatrix[,S+1] + AllSpAbunds[,s]
      }else{
     SpMatrix[,i] <- AllSpAbunds[,s]
     i <- i + 1}
    }else{next}
  }
   
   SpNames <- c(AllSpNames[SpToKeep],"rare")
   names(SpMatrix) <-  SpNames
   Intra <- ifelse(SpNames == "conspecific", 1, 0)
   SpMatrix[,which(SpNames ==specific.focal)] <- SpMatrix[,which(SpNames == specific.focal)] - SpMatrix[,"conspecific"]  # remove conspecific obs from focal family
  
   SpMatrix <- as.matrix(SpMatrix)
  # max fecundity
  Nmedian <- c(SpMatrix[which.max(Fecundity),])

  # Upper bound intrinsic fecundity
 Fmean <- ceiling(log(mean(Fecundity)))
  
  
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
                  S=S+1,
                  Fmean =Fmean  ,
                  Nmedian = Nmedian ,
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
  list.init <- function(...)list(lambda= array(as.numeric(Fmean ), dim = 1))
                                 #N_opt = array(as.numeric(Nmax), dim = S+1))
  
  

  if( run.stan == 1){                               
  FinalFit <- stan(file = "code/DensityFunct_Final.stan", 
                   data = DataVec,
                   init=  list.init,
                   #init =  "random",
                   warmup= 500,
                   iter = 1000, 
                   init_r = 2,
                   chains = 4,
                   control=list(max_treedepth=15),
                   seed= 1644)
  
  
  save(file= paste0("results/stan/FinalFit_",grouping,"_",Code.focal,"_",alpha.function,".rds"),
       FinalFit)
  
  if(length(rstan::extract(FinalFit)) == 0){ # restrict less niital values
   FinalFit <- stan(file = "code/DensityFunct_Final.stan", 
                     data = DataVec,
                     init =  "random",
                     warmup= 500,
                     iter = 1000, 
                     init_r = 2,
                     chains = 4,
                     control=list(max_treedepth=15),
                     seed= 1644)
   save(file= paste0("results/stan/FinalFit_",grouping,"_",Code.focal,"_",alpha.function,".rds"),
        FinalFit)
    
  }
  }
  #load(paste0("results/stan/FinalFit_",Code.focal,"_",alpha.function,".rds"))
  
 load(paste0("results/stan/FinalFit_",grouping,"_",Code.focal,"_",alpha.function,".rds"))

 FinalPosteriors <- rstan::extract(FinalFit)
  
  print("Final Fit done")
  
  #---- 3.3. Final fit posterior check and behavior checks---- 
  if( run.diagnostic == 1){  
    if(function.int == 1){par <- c('lambda',
                                   'alpha_initial')
    }
    if(function.int ==2){par <- c('lambda',
                                  'alpha_initial',
                                  'alpha_slope')
    }
    if(function.int >2){par <- c('lambda',
                                 'c','alpha_initial',
                                 'alpha_slope','N_opt')
    }
    
  ##### Diagnostic plots and post prediction 
  pdf(paste0("figures/stan/FinalFit_",grouping,"_",Code.focal,"_",alpha.function,".pdf"))
  # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
  source("code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
  # check the distribution of Rhats and effective sample sizes 
  ##### Posterior check
  stan_post_pred_check(FinalPosteriors,"F_hat",Fecundity,
                       paste0("results/stan/PostFec_",grouping,"_",Code.focal,"_",alpha.function,".csv.gz")) 
  
  #log_lik_2 <- loo::extract_log_lik(FinalFit, 
  #                                  parameter_name = "F_sim", 
  #                                  merge_chains = F)
  
  #r_eff <- loo::relative_eff(exp(log_lik_2), cores = 2) 
  
 # loo_1 <- loo::loo(log_lik_2 , r_eff = r_eff, cores = 2)
  
 # print(loo_1)
  
  # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
  hist(summary(FinalFit)$summary[,"Rhat"],
       main = paste("Finat Fit: Histogram of Rhat for",
                    Code.focal," and ",alpha.function))
  hist(summary(FinalFit)$summary[,"n_eff"],
       main = paste("Finat Fit: Histogram of Neff for",
                    Code.focal," and ",alpha.function))
  
  # plot the corresponding graphs
  trace <- stan_trace(FinalFit, pars=par,
                      inc_warmup = TRUE)
  print(trace)
  dens <- stan_dens(FinalFit, 
                    pars=par)
  print(dens)
  splot <- stan_plot(FinalFit, 
                     pars=par)
  print(splot)
  sampler_params <- get_sampler_params(FinalFit, inc_warmup = TRUE)
  summary(do.call(rbind, sampler_params), digits = 2)
  pairs(FinalFit, pars = par)
  
  # Next check the correlation among key model parameters and identify any
  #pairs(FinalFit, pars = c("lambda",'alpha_initial','alpha_slope','c'))
  
  
  dev.off()
  }
  #---- 3.3. Extract coefficients ----
  
  assign(paste0("Parameters_",Code.focal,"_",alpha.function),
       list(DataVec = DataVec,
            SpNames=SpNames,
       alpha_value= FinalPosteriors$alpha_value,      
       lambda =  FinalPosteriors$lambda_ei,
       alpha_slope = FinalPosteriors$alpha_slope,
       alpha_init =FinalPosteriors$alpha_init,
       alpha_c = FinalPosteriors$c,
       N_opt = FinalPosteriors$N_opt
       ))
  
  save(list =paste0("Parameters_",Code.focal,"_",alpha.function),
       file = paste0("results/stan/Parameters_",grouping,"_",Code.focal,"_",alpha.function,".RData"))
  
  #---- 3.3. Extraction fecundity---
  
  assign(paste0("Fsim_",Code.focal,"_",alpha.function),
         list(DataVec = DataVec,
              F_sim=FinalPosteriors$F_sim
         ))
  
  save(list =paste0("Fsim_",Code.focal,"_",alpha.function),
       file = paste0("results/stan/Fsim_",grouping,"_",Code.focal,"_",alpha.function,".RData"))

}
}

#---- 4. Figures ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---- 4.0. check models ----
ModelCheck <- NULL
grouping<- "family"
for(Code.focal in "LARO"){#focal.levels
    for (function.int in c(1:4)){ # c(1:4)
      print(paste(Code.focal,", function",function.int))
      
      alpha.function <- paste0("function_",function.int)
      
      load(paste0("results/stan/FinalFit_",grouping,"_",Code.focal,"_",alpha.function,".rds"))
      
      mc <- data.frame(focal =Code.focal, 
                        function.int = function.int,
                        Rhat = max(summary(FinalFit)$summary[,"Rhat"],na.rm =T),
                        Neff = min(summary(FinalFit)$summary[,"n_eff"],na.rm = T))
      ModelCheck <- bind_rows(ModelCheck,mc)
      remove(FinalFit)
    }
  }
View(ModelCheck)

write.csv(ModelCheck,
          paste0("results/ModelCheck_LARO_",grouping,".csv"))

#---- 4.1 Alpha distribution ----

df_alpha_all <- data.frame(focal=NA)
for(Code.focal in "LARO" ){ #focal.levels
  for (function.int in c(1:4)){
    df_alpha <- data.frame()
    function.vec <- c(0,0,0,0)
    function.vec[function.int] <- 1
    alpha.function <- paste0("function_",which(function.vec==1))

    load(paste0("results/stan/Parameters_",grouping,"_",Code.focal,"_",alpha.function,".RData"))
    parameters <- get(paste0("Parameters_",Code.focal,"_",alpha.function))
    neighbours.vec <- c(parameters[["SpNames"]])
    
    N <- parameters[["DataVec"]]$N
    Nmax <- parameters[["DataVec"]]$Nmax
    df_alpha_init <- data.frame(observation = c(1:N),
                                 focal = Code.focal,
                                 function.int = function.int)
     df_alpha_init <- NULL
     for(n in 1:length(neighbours.vec)){
       df_alpha_init_n <- data.frame(focal = Code.focal,
                                     Nmax=parameters[["N_opt"]][,n],
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

levels(as.factor(df_alpha_all$neigh))
source("code/1.1.PopProjection_toolbox.R")

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
          paste0("results/df_funct_alpha_LARO_",grouping,".csv.gz"))

df_funct_alpha <- data.table::fread("results/df_funct_alpha_LARO_family.csv.gz")

cbp2 <- c("#000000", "#E69F00", "#56B4E9","#F0E442","#009E73",
          "#CC79A7", "#0072B2", "#D55E00")

cb2_focal <- c("#F0E442","#009E73",
               "#CC79A7", "#0072B2")
funct.col <- c("#AA4499", "#DDCC77","#88CCEE", "#44AA99")
plot.list.alpha <- list()


# percentage of positive interaction with function 4
100*nrow(df_funct_alpha[which(df_funct_alpha$function.int == 4 &
                       df_funct_alpha$alpha_value > 0),])/
  nrow(df_funct_alpha[which(df_funct_alpha$function.int == 4),])


# Only LARO
df_funct_alpha <- df_funct_alpha %>%
  mutate(function.name = case_when(function.int==1 ~"1.Traditional",
                                   function.int==2 ~"2.Linear",
                                   function.int==3 ~"3.Exp",
                                   function.int==4 ~"4.Sigmoid"))


dummy <- data.frame(uplimit=c(0.060,0.5,1.6,1),
                    downlimit=c(-0.1,-1,-0.5,-0.5),
                    function.name = c("1.Traditional","2.Linear","3.Exp","4.Sigmoid"),
                    stringsAsFactors=FALSE)

plot_LAROalpha_list <- list()
family.neigh <- levels(as.factor(df_funct_alpha$neigh))
family.neigh <- family.neigh[!family.neigh=="conspecific"]
for( i in c("1.Traditional","2.Linear","3.Exp","4.Sigmoid")){
  df <- df_funct_alpha %>%
    dplyr::filter(focal == "LARO" & function.name==i)%>%
    mutate(neigh= factor(neigh, 
                            levels = c(family.neigh,"conspecific")))
  plot_LAROalpha_list[[i]] <- ggplot(df,
                                       aes(x=density, y= alpha_value,
      color=neigh,fill=neigh,size=neigh,alpha=neigh)) +
    geom_line(stat="smooth",method = 'gam',
              se = TRUE,level =0.95,
              alpha=0.7,
              size=3) +
    #stat_smooth(method = 'gam',se = TRUE,level =0.95,aes(size=neigh,alpha=neigh)) +
    geom_hline(yintercept=0, color="black", linetype="dashed",
               size=2) +
    scale_color_manual("Neighbours identity",values=rev(cbp2)) + 
    scale_fill_manual("Neighbours identity",
                      values=rev(cbp2)) + 
    scale_size_manual("Neighbours identity",
                      values=c(3, rep(1.5,length(family.neigh)))) + 
    scale_alpha_manual("Neighbours identity",
                       values=c(1, rep(0.8,length(family.neigh)))) + 
    theme_bw() +
    labs(title = i ,
      y= "Per capita effect of neighbours on focal",
      x="density of neighbours") +
    rremove("ylab") + rremove("xlab") + 
    scale_x_continuous( expand= c(0,0),
                        minor_breaks = NULL,
                        breaks=c(0,5,10),
                        labels=c("0","5","10")) +
    scale_y_continuous(expand= c(0,0),
                       minor_breaks = NULL,
                       breaks=c(dummy$downlimit[which(dummy$function.name==i)],
                                0,
                                dummy$uplimit[which(dummy$function.name==i)])) +
    guides(color= "none",fill="none",size="none",alpha="none") +
    coord_cartesian( xlim = NULL, 
                     #ylim=c(-1,1),
                     ylim = c(dummy$downlimit[which(dummy$function.name==i)],
                     dummy$uplimit[which(dummy$function.name==i)]),
                     expand = TRUE, default = FALSE, clip = "on") +
    theme(plot.title = element_blank(),#element_text(color="#000000",
                                    #vjust = 0,
                                    #size=24),
          
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text =element_text(color="black",size=20),
          legend.key= element_rect(fill = "white"),
          legend.text=element_text(size=20),
          legend.title=element_blank(),
          axis.text.x = element_text(color="black",size=16),
          axis.text.y = element_text(color="black",
                                     size=16),
          legend.key.size = unit(1, 'cm'),
          plot.margin = unit(c(1,1,0,0), 'lines')) 
  
}

plot_LAROalpha <- ggarrange( plotlist = plot_LAROalpha_list,
                              ncol=1,align = ("v"))

plot_LAROalpha <- annotate_figure(plot_LAROalpha , 
                                  left = textGrob("Per capita effect of neighbours on focal", 
                                                  rot = 90, vjust = 0.5, #hjust=1, 
                                                  gp = gpar(fontsize=20,cex = 1.3)),
                                  bottom = textGrob("Density of neighbours",  
                                                    gp = gpar(fontsize=20,cex = 1.3))) 
plot_LAROalpha
ggsave(plot_LAROalpha,
       file = "figures/NatData_plot_LARO.pdf")

#---- 4.2. Population Projections ----

df_param_all <- data.frame(focal=NA)

for(Code.focal in "LARO"){
  for (function.int in c(1:4)){
    #if (Code.focal == "TROR") next
    df_param <- data.frame()
    function.vec <- c(0,0,0,0)
    function.vec[function.int] <- 1
    alpha.function <- paste0("function_",which(function.vec==1))
    
    load(paste0("results/stan/Parameters_family_",Code.focal,"_",alpha.function,".RData"))
    parameters <- get(paste0("Parameters_",Code.focal,"_",alpha.function))
    N <- parameters[["DataVec"]]$N
    Nmax <- parameters[["DataVec"]]$Nmax
    neighbours.vec <- c(parameters[["SpNames"]])
    df_param_init <- data.frame(observation = c(1:N),
                                focal = Code.focal,
                                function.int = function.int)
    df_param_init <- NULL
    for(n in 1:length(neighbours.vec)){
      df_param_init_n <- data.frame(focal = Code.focal,
                                    Nmax=parameters[["N_opt"]][,n],
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

community_id_df <- data.table::fread("/Users/lisabuche/Documents/Projects/Perenjori/results/community_id_df_short.csv")
str(community_id_df)
levels(as.factor(community_id_df$family))

community_id_df_LARO <- community_id_df %>% 
  dplyr::filter(code=="LARO" ) %>%
  aggregate(count ~ code + year + id.plot + collector +  family + scale.weight, sum) %>%
  mutate(count=count/16)

head(community_id_df_LARO)
ggplot(community_id_df_LARO, aes(x=count)) +geom_density()

fam.to.keep <- levels(as.factor(df_param_all$neigh))
fam.to.keep <- fam.to.keep[which(!fam.to.keep =="rare")]

community_id_df_fam <- community_id_df %>% 
  mutate(groups = case_when(!family %in% fam.to.keep ~ "rare",
                           T ~ family)) %>%
  #dplyr::filter(family %in% levels(as.factor(df_param_all$neigh))) %>%
  aggregate(count ~ id.plot + collector+ year + groups + scale.weight, sum) %>%
  mutate(count=count/16) %>%
  spread(groups,count)
str(community_id_df_fam)

df_projection <- NULL
gens = 20
for(Code.focal in "LARO"){
  for(i in 1:100){
  state_df <- data.frame("araliaceae"= round(sample(community_id_df_fam$araliaceae[!is.na(community_id_df_fam$araliaceae)], gens)),
                         "asteraceae"=round(sample(community_id_df_fam$asteraceae[!is.na(community_id_df_fam$asteraceae)], gens)),
                         "crassulaceae"=round(sample(community_id_df_fam$crassulaceae[!is.na(community_id_df_fam$crassulaceae)], gens)),
                         "goodeniaceae"=round(sample(community_id_df_fam$goodeniaceae[!is.na(community_id_df_fam$goodeniaceae)], gens)),
                         "montiaceae"=round(sample(community_id_df_fam$montiaceae[!is.na(community_id_df_fam$montiaceae)], gens)),
                         "poaceae"=round(sample(community_id_df_fam$poaceae[!is.na(community_id_df_fam$poaceae)],gens)),
                         "conspecific" = round(sample(community_id_df_LARO$count[which(community_id_df_LARO$year == 2022)],1,replace = T)),
                         "rare"=round(sample(community_id_df_fam$rare[!is.na(community_id_df_fam$rare)],gens)))
  for (function.int in c(1:4)){
    param.df <- df_param_all[which(df_param_all$function.int==function.int & 
                                     df_param_all$focal==Code.focal ),]
    param.df$g = 0.7
    param.df$s = 0.9
   df_projection_n <- Ricker_solution_NatData(gens=gens,
                                             state=state_df,
                                             pars = param.df)

  df_projection_n$sim <- i 
  df_projection_n$function.int <- function.int
  df_projection_n$focal <- Code.focal
  df_projection_n$year <- c(2023:2042)
  df_projection <- bind_rows(df_projection, df_projection_n)
  }
  }
 }

df_projection <- df_projection %>%
  mutate(function.name = case_when(function.int==1 ~"1.Traditional",
                                   function.int==2 ~"2.Linear",
                                   function.int==3 ~"3.Exp",
                                   function.int==4 ~"4.Sigmoid")) %>%
  gather(all_of(levels(as.factor(df_param_all$neigh))), 
         key="species", value="abundance")   

community_id_df_projectiongraph  <- community_id_df %>% 
  dplyr::filter(genus =="Lawrencella"|family %in% levels(as.factor(df_param_all$neigh))) %>%
  mutate(family = case_when(genus =="Lawrencella" ~"conspecific",
                            T~family)) %>%
  aggregate(count ~ id.plot + collector +  family + year + scale.weight, sum) %>%
  mutate(count=count/16) %>% 
  rename("species" = family,
        "abundance"= count) %>%
  select(year, species, abundance)
  
community_id_df_projectiongraph.1 <- community_id_df_projectiongraph %>%
  mutate(function.name ="1.Traditional")
community_id_df_projectiongraph.2 <- community_id_df_projectiongraph %>%
  mutate(function.name ="2.Linear") %>%
  bind_rows(community_id_df_projectiongraph.1)
community_id_df_projectiongraph.3 <- community_id_df_projectiongraph %>%
  mutate(function.name ="3.Exp") %>%
  bind_rows(community_id_df_projectiongraph.2)
community_id_df_projectiongraph.4 <- community_id_df_projectiongraph %>%
  mutate(function.name ="4.Sigmoid")%>%
  bind_rows(community_id_df_projectiongraph.3)

allyear_communityprojection <- bind_rows(df_projection,
                                         community_id_df_projectiongraph.4 )


dummy <- data.frame(uplimit=c(80,60,40,20),
                    ylab=c(70,52,40,16),
                    significance=c("*","","**","***"),
                    function.name = c("1.Traditional","2.Linear","3.Exp","4.Sigmoid"),
                    stringsAsFactors=FALSE)

family.neigh <-levels(as.factor(allyear_communityprojection$species))
family.neigh <-  family.neigh[!family.neigh =="conspecific"]

plot_projection_list <- list()
for( i in c("1.Traditional","2.Linear","4.Sigmoid")){ #"3.Exp"
  df <- allyear_communityprojection %>%
    dplyr::filter(function.name==i) %>%
    mutate(species = factor(species, 
                            levels = c(family.neigh,"conspecific")))
  plot_projection_list[[i]] <- ggplot(df) +
    annotate("rect", xmin=2010,xmax=2022,
            ymin=-5,ymax=1000,
             fill="lightgrey",alpha=0.3) + 
    annotate("text", x=2020, y=dummy[which(dummy$function.name==i),"ylab"],
               label=dummy[which(dummy$function.name==i),"significance"],
               color="black",size=16) + 
    stat_summary(data=df[which(df$species=="conspecific"),],
                 aes(y=abundance, x=year),
                 color="black",
                 alpha=0.9,size=1.5,
                 fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange") +
    stat_summary(aes(y=abundance, x=year,
                     group=species,color=species,
                     size =species, alpha=species),
                 fun.y = mean,
                 geom = "point") +
    stat_summary(aes(y=abundance, x=year,
                     group=species,color=species,
                     size =species, alpha=species),
                 fun.y = mean,
                 geom = "line") +
   #scale_color_manual("",values=rev(cbp2),
    #                  labels=c("Araliaceae","Asteraceae","Crassulaceae",
    #                           "Goddeniaceae","Montiaceae","Poaceae","Rare families",
    #                           "LARO")) + 
   scale_color_manual("",values=rev(c("black",rep("grey",
                                                  times=length(family.neigh))))) +#values=cbp2) + 
   scale_size_manual("",values=rev(c(1.5, rep(1.2,length(family.neigh)))),
                     labels=c("Araliaceae","Asteraceae","Crassulaceae",
                              "Goddeniaceae","Montiaceae","Poaceae","Rare families",
                              "LARO")) + 
   scale_alpha_manual("",values=rev(c(0.9, rep(0.6,length(family.neigh)))),
                      labels=c("Araliaceae","Asteraceae","Crassulaceae",
                               "Goddeniaceae","Montiaceae","Poaceae","Rare families",
                               "LARO")) + 
    theme_bw() +
    labs(title = i,
         y= "",
         x="") +
    rremove("ylab") + rremove("xlab") + 
    guides(color= guide_legend("",position="bottom",
                               direction="horizontal",
                          byrow = TRUE,
                          nrow = 2)) +
    coord_cartesian( xlim = NULL, ylim = c(-1,dummy$uplimit[which(dummy$function.name==i)])) + #ylim =c(-5,80)) + #
    scale_x_continuous(expand= c(0,0),minor_breaks = NULL) +
    scale_y_continuous(expand= c(0,0),
                       minor_breaks = NULL,
                       breaks=c(0,20,40,60,80),
                       labels=c(0,20,40,60,80)) +
    theme_bw() +
      facet_wrap(.~function.name, strip.position="right") +
    theme(plot.title = element_blank(), #element_text(color="#000000",
                                    #vjust = 0,size=24),
          strip.background = element_blank(),
          strip.placement = "right",
          strip.text =element_text(color="black",
                                   size=22),
          legend.key= element_rect(fill = "white"),
          legend.text=element_text(size=16),
          legend.title=element_blank(),
          axis.text.x = element_text(color="black",size=20),
          axis.text.y = element_text(color="black",size=20),
          legend.key.size = unit(1, 'cm'),
          plot.margin = unit(c(1,1,0,0), 'lines'))
  
}
library(ggpubr)
legend.plot <- as_ggplot(ggpubr::get_legend(  plot_projection_list[[i]]))

library(grid)
plot_projection <- ggpubr::ggarrange( plotlist = plot_projection_list,
           ncol=1,
           align = c("v"),
           common.legend = F,
           legend="none")


plot_projection  <- annotate_figure(plot_projection  , 
                              left = textGrob("Density at 25 x 25cm scale", 
                                              rot = 90, vjust = 0.5,
                                              gp = gpar(fontsize=20,cex = 1.3)),
                              #bottom = textGrob("Time", 
                              #                  #vjust = -0.5,
                              #                  gp = gpar(fontsize=20,cex = 1.3)),
                              top= textGrob("Past observations                    Predicted density", 
                                            rot = 0, hjust = 0.55,
                                            gp = gpar(fontsize=20,cex = 1.3))
                              ) +
  theme(plot.margin = unit(c(0,0,0,2), 'lines'))

plot_projection 

plot_projection_alpha <- ggarrange(plot_LAROalpha,
                                   plot_projection,
                                   ncol=2, 
                                   align = c("h"),
                                   widths=c(1,2),
          legend = "none", common.legend = F) 

plot_projection_alpha <- ggarrange(plot_projection_alpha ,
                                   legend.plot,
                                   ncol=1,nrow=2,
                                   align = c("v"),
                                   heights=c(8,1),
                                   legend = "none", 
                                   common.legend = F) 
plot_projection_alpha 
ggsave(plot_projection_alpha , 
       file = "figures/NatData_plot_projection.pdf")


plot_projection

ggsave(plot_projection , 
       file = "figures/NatData_plot_poster.pdf")

#---- 4.3. Compare model for each focal  ----
# reference : Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing. 27(5), 1413–1432. :10.1007/s11222-016-9696-4. online, arXiv preprint arXiv:1507.04544.
model.loo <- list()

for(Code.focal in "LARO"){
  for (function.int in 1:4){ # c(1:4)
    if(ModelCheck[which(ModelCheck$function.int==function.int &
                       ModelCheck$focal==Code.focal),"Rhat"] < 0 | 
      ModelCheck[which(ModelCheck$function.int==function.int &
                         ModelCheck$focal==Code.focal),"Neff"] < 100) next
      
      print(paste(Code.focal,", function",function.int))
      
      load(paste0("results/stan/FinalFit_family_",Code.focal,"_function_",function.int,".rds"))
      # Extract pointwise log-likelihood
      # using merge_chains=FALSE returns an array, which is easier to 
      # use with relative_eff()
      log_lik <- loo::extract_log_lik(FinalFit, 
                                      parameter_name = "F_sim", 
                                      merge_chains = F)
      #as of loo v2.0.0 we can optionally provide relative effective sample sizes
      # when calling loo, which allows for better estimates of the PSIS effective
      # sample sizes and Monte Carlo error
      r_eff <- loo::relative_eff(log_lik, cores = 3) 
      # preferably use more than 2 cores (as many cores as possible)
      # will use value of 'mc.cores' option if cores is not specified
      model.loo[[paste0(Code.focal,"_function_",function.int)]] <-  loo::loo(log_lik,
                                                                            threshold=0.7,
                                                                            r_eff = r_eff, 
                                                                            cores = 3)
      remove(FinalFit)
  }
}

Se_loo_model <- NULL
for(Code.focal in "LARO"){
  Se_loo_model <- loo::loo_compare(model.loo[[paste0(Code.focal,"_function_1")]], 
                        model.loo[[paste0(Code.focal,"_","function_2")]],
                        model.loo[[paste0(Code.focal,"_","function_3")]],
                        model.loo[[paste0(Code.focal,"_","function_4")]])
  

  Se_loo_model <-as.data.frame(print(Se_loo_model, simplify = FALSE, digits = 3))
  
  }
Se_loo_model<- Se_loo_model %>%
  rownames_to_column(var="model") #%>%
  #separate(col=function.int, 
  #         into=c("model","function.int"),sep="l") %>%
  #dplyr::subset(-model)


write.csv(Se_loo_model,
          paste0("results/LARO_Se_loo_model.csv"))

#---- 4.3. Figure of fecundity estimation distribution ----
source("code/stan_modelcheck_rem.R")
function.names <- c("1.Traditional","2.Linear","3.Exp","4.Sigmoid")
  pdf(paste0("figures/NatData_PostFecundity_distribution.pdf")) 
  par(oma=c(0,0,0,0), mar=c(2.25,2.5,2,1))
  layout(mat = matrix(
    1:4,
    byrow=TRUE,
    nrow = 1,
    ncol = 4
  ),
  heights = rep(1,8),
  widths = rep(4,8)
  )
  
  for(Code.focal in "LARO"){ #,"j"
    for (function.int in c(1:4)){ # c(1:4)
     #col2 <- c("black","#CC79A7","#E69F00","#009E73") # function 
      col2 <- wes_palette("FantasticFox1", n = 5)
      #col1 <- c("blue","red") # focal species
      # loo values
      value.se <- round(Se_loo_model$looic[which(Se_loo_model$model== paste0("model",function.int))],
                        digits=2)
      # fecundity obs
      load(paste0("results/stan/Fsim_family_",Code.focal,"_function_",function.int,".RData"))
      # fecundity generated
      Fec_df <- data.frame(Obs = 1:get(paste0("Fsim_",Code.focal,"_function_",function.int))[["DataVec"]]$N,
                           Fec =get(paste0("Fsim_",Code.focal,"_function_",function.int))[["DataVec"]]$Fecundity)
      
      
      load(paste0("results/stan/FinalFit_family_",Code.focal,"_function_",function.int,".rds"))
      FinalPosteriors <- rstan::extract(FinalFit)
      stan_post_pred_check_all(FinalPosteriors,"F_hat",
                               Fec_df$Fec,
                               function.names[function.int],
                               "black",
                               col2[1+function.int],value.se
                               ) 
      
      
    }
  }
  dev.off()




