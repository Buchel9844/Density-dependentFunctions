#### coexistence evaluation for parameter drawn from entire range of potemntial value
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
library(odeintr)
library(brms)
library(bbmle)
library(reshape2)
library(ggthemes)
library(tsvr) # library for the timescale specific variance ratio

library(ppcor) # partial correlation value

 # for ##########################################################################################################
# 1. Compute Low- Density Growth Rates for wide range of parameters
##########################################################################################################
source("code/PopProjection_toolbox.R")

# define parameters
nsims <- 500
species <- 2
function.int <- c(1:4)

# run sim
set.seed(1608)

# Set parameters lists
alpha.coexistence <- function(species){
  a.vec <- sort(rtruncnorm(n=species*species, a=-1, b=0, mean=0, sd=0.2))
  a.vec <- c(a.vec[1],a.vec[3],a.vec[4],a.vec[2])
  mat <- matrix(ncol=species, nrow = species, # a_slope, only negative values
         data = a.vec,
         dimnames = list(c("i", "j"),
                         c("i", "j")))
  return(mat)
}
params <-list()
for( i in 1:nsims){
    params[[i]] <- list(
                # stable parameters
                    sim= i, 
                    g = c(gi=1,gj=1), # germination rate of seeds
                    s = c(si= 1,sj= 1), # survival rate of seeds
                    lambda = c(lambdai= 1 ,lambdaj= 1), # intrinsic fecundity
                    Nmax =matrix(ncol=species, nrow = species, # N optimal, only positive values
                                 data = round(abs(runif(species*species, max=10, min = 0))),
                                 dimnames = list(c("i", "j"),
                                                 c("i", "j"))),
                # variable parameters
                    a_slope = matrix(ncol=species, nrow = species, # a_slope, only negative values
                                  data = runif(n=species*species, min=-1, max=0),
                                  dimnames = list(c("i", "j"),
                                                  c("i", "j"))),
                    
                    a_initial = alpha.coexistence(species),
              
                    c =matrix(ncol=species, nrow = species, # a_slope, only positive values
                                     data = rtruncnorm(n=species*species, a=0, b=1, mean=0, sd=1),
                              dimnames = list(c("i", "j"),
                                              c("i", "j")))
                    )
    if(params[[i]]$a_initial[1,1]*params[[i]]$a_initial[2,2] < 
       params[[i]]$a_initial[1,2]*params[[i]]$a_initial[2,1]){
      nsims = nsims - 1}

}


df.sim  <- NULL

t.num = 100 # number of generation
for(i in 1:nsims){
  for( function.int in 1:4){
    print(paste0("int ", i,"for funct ",function.int))
    df.n <- NULL
    p <- params[[i]]
    df.n <- as.data.frame(params[[i]]) %>%
      mutate(focal = c("i","j"))  %>%
      melt(id.vars="focal")   %>%
      mutate( vars = paste(variable ,focal,sep=".")) %>%
      dplyr::select( value, vars) %>%
      spread(key=vars, value=value) %>%
      mutate(function.int  = function.int ) # attention!!! identity of the focal is the last suffixes
    
    df.inv <- NULL
    df.inv <- GrowthSimInv(par.dat = p, t.num  = t.num,
         function.int) %>%# simulation of invasion for both species
      bind_cols(df.n) 
  
    df.initc <- Ricker_solution_ODE(state = c(1,1), pars= p, gens=t.num*10,
                        function.int) %>% # simulation of initial condition = 1,1
      mutate(invader = "both") %>%
      bind_cols(df.n) %>% 
      bind_rows(df.inv)

    df.sim <- bind_rows(  df.sim, df.initc)

  }
}

# save .csv

write_csv(x = df.sim , col_names = TRUE, 
          file = paste0("results/df.sim.csv.gz"))


##########################################################################################################
# 2. Visualise population dynamics
##########################################################################################################
sim = 2
t.num= 100
#---- Abundance through time for all scenario----
example.abundance.scenarios <- df.sim[which(df.sim$sim.i == sim &
               df.sim$function.int == 4),] %>%
  gather(Ni, Nj, key=species, value=abundance) %>%
  ggplot(aes(y=abundance, x=time)) + #geom_point() +
  geom_smooth(aes(y=abundance,color=species,fill=species),
              alpha=0.2,size=0.5, linetype="dashed") +
  geom_line(aes(y= abundance,color=species), alpha=0.8) +
  xlim(c(0,100)) +
  theme_bw() + facet_wrap(as.factor(invader)~., scales="free") + 
  scale_colour_colorblind() +
  scale_fill_colorblind() +
  labs(title="abundances Ni and Nj over time")

ggsave(example.abundance.scenarios,
       file = "figures/example.abundance.scenarios.pdf")

#---- Abundance through time----
example.oscillatory.state <- df.sim[which(df.sim$sim.i == sim &
               df.sim$invader == "both" &
                     (df.sim$function.int == 4 | 
                        df.sim$function.int == 1)),] %>%
  gather(Ni, Nj, key=species, value=abundance) %>%
  ggplot(aes(y=abundance, x=time)) + #geom_point() +
  geom_smooth(aes(y=abundance,color=species,fill=species),alpha=0.2,size=0.5, linetype="dashed") +
  geom_line(aes(y= abundance,color=species),alpha=0.8) +
  theme_bw() + facet_wrap(function.int~.) + 
  scale_colour_colorblind() +
  xlim(c(0,100)) +
  scale_fill_colorblind() +
  labs(title="abundances Ni and Nj over time when 
       Ni = Nj = 0 at t=0")

ggsave(example.oscillatory.state,
       file = "figures/example.oscillatory.state.pdf")

#---- GR through time----

GRWR.list <- list()

for(function.int in c(1:4)){
  par.dat <- params[[sim]]

GRWR.list[[function.int]] <- df.sim[which(df.sim$sim.i == 20 &
                                            df.sim$invader == "both" &
                                            df.sim$function.int == function.int),] %>%
  gather(dNi, dNj, key=species, value=GR) %>%
  ggplot(aes(y=GR, x=time))  +
  xlim(c(0,100)) + 
  geom_line(aes(y=GR,color=species),alpha=0.7,size=1)  + 
  geom_vline(aes(xintercept=1),linetype="dashed",color="black",alpha=0.7) +
  theme_bw() + scale_colour_colorblind() +
  guides(alpha="none") +
  labs(title=paste0("GRWR at each time step for Ni and Nj \n for sim ",sim," and function ",function.int)) 
}
example.GR.all.timesteps <- ggarrange(plotlist = GRWR.list, ncol=2, nrow=2,common.legend = T,
          legend = "right")

ggsave(example.GR.all.timesteps,
       file = "figures/example.GR.all.timesteps.pdf")

#---- Ni dependency on Nj ----
NINJ.list <- list()

for(function.int in c(1:4)){
  par.dat <- params[[sim]]
  
  NINJ.list[[function.int]] <- df.sim[which(df.sim$sim.i == sim &
               df.sim$invader == "both" &
               df.sim$function.int == function.int ),]  %>%
  ggplot(aes(Ni, Nj)) + 
  geom_path(aes(colour = as.numeric(time))) +
  facet_wrap(function.int~.) + 
  theme_bw() +
 labs(title=paste0("Abundance at each time step 
 for Ni and Nj \n 
 for sim ",sim,
                   " and function ",function.int)) 
  

}
example.oscillatory.state.of.abundances <- ggarrange(plotlist =   NINJ.list, ncol=2, nrow=2,common.legend = T,
          legend = "right")


ggsave(example.oscillatory.state.of.abundances,
       file = "figures/example.oscillatory.state.of.abundances.pdf")

#---- GR of Ni dependency on GR of Nj ----

example.oscillatory.state.of.GR <- df.sim[which(df.sim$sim.i == sim &
               df.sim$invader == "both" &
               df.sim$function.int == function.int),]  %>%
  ggplot(aes(dNi, dNj)) + 
  geom_path(aes(colour = as.numeric(time))) +
  facet_wrap(function.int~.) + 
  theme_bw() 

ggsave(example.oscillatory.state.of.GR,
       file = "figures/example.oscillatory.state.of.GR.pdf")

example.oscillatory.boundaries <- ggarrange( ncol=2, nrow=1,
df.sim[which(df.sim$sim.i == sim &
               df.sim$invader == "both"),]  %>%
  mutate(Nj2 = c(Nj[-1],NA))%>%
  ggplot(aes(Nj, Nj2)) + 
  geom_path(aes(colour = as.numeric(time))) +
  facet_wrap(function.int~., scales="free") + 
  labs(title="Nj boundaries", x= "N(t+1)", y = " N(t)") +
  theme_bw(),
df.sim[which(df.sim$sim.i == sim &
               df.sim$invader == "both"),]  %>%
  mutate(Ni2 = c(Ni[-1],NA)) %>%
  ggplot(aes(Ni, Ni2)) + 
  geom_path(aes(colour = as.numeric(time))) +
  facet_wrap(function.int~., scales="free") +
  labs(title="Ni boundaries", x= "N(t+1)", y = " N(t)") +
  theme_bw() 
)

ggsave(example.oscillatory.boundaries,
       file = "figures/example.oscillatory.boundaries.pdf")

