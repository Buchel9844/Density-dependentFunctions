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
library(data.table) # write csv.gz
library(ppcor) # partial correlation value

###########################################################################################################
# 1. # 1. Compute Low- Density Growth Rates for wide range of parameters without coexistence restriction
##########################################################################################################
source("code/PopProjection_toolbox.R")

# define parameters
nsims <- 500
species <- 2
function.int <- c(1:4)
t.num = 100
gen = seq(from=0, to=t.num , by=1)
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
white.noise <- function(t.num){
  e_noise = rnorm(t.num+1,mean=0, sd=0.1) 
  for( n in 1:sample(2:10)[1]){
    e_noise[n] <- rnorm(1,mean=0,
                        sd=abs(rnorm(1, mean=0, sd=1)))
  }
  return( e_noise)
}

params <-list()
i=1
for( n in 1:10000){
  if(i==nsims+1) next
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
                    
                    a_initial = matrix(ncol=species, nrow = species, # a_slope, only negative values
                                       data = rtruncnorm(n=species*species, mean=-0.2, sd=0.5, a=-1, b=1),
                                       dimnames = list(c("i", "j"),
                                                       c("i", "j"))),
              
                    c =matrix(ncol=species, nrow = species, 
                                     data =runif(n=species*species, min=-1, max=0),
                              dimnames = list(c("i", "j"),
                                              c("i", "j"))),
              
                    e_seasonal = sin((2*pi/20)*gen)*abs(rnorm(1,mean=0, sd=0.1)),
                    e_noise = white.noise(t.num)
                  
                    )
    if(params[[i]]$a_initial[1,1] > 0 | params[[i]]$a_initial[2,2] > 0){
      i = i}else{
        i = i + 1
      }
}


df.sim  <- NULL

t.num = 100 # number of generation
set.seed(1608)
for(i in 1:nsims){
  for( function.int in 1:4){
    for(add_external_factor in c("none","season","noise")){
    print(paste0("int ", i,"for funct ",function.int, add_external_factor))
    df.n <- NULL
    p <- params[[i]]
    df.n <- as.data.frame(params[[i]][c("sim","Nmax","a_slope","a_initial","c")]) %>%
      mutate(focal = c("i","j"))  %>%
      reshape2::melt(id.vars="focal")   %>%
      mutate( vars = paste(variable ,focal,sep=".")) %>%
      dplyr::select( value, vars) %>%
      spread(key=vars, value=value) %>%
      mutate(function.int  = function.int ) # attention!!! identity of the focal is the last suffixes
    
    df.inv <- NULL
    df.inv <- GrowthSimInv(par.dat = p, t.num  = t.num,
         function.int,
         add_external_factor) %>%# simulation of invasion for both species
      bind_cols(df.n) 
   
    df.initc <- Ricker_solution(state = c(5,5), pars= p, gens=t.num,
                        function.int,add_external_factor) %>% # simulation of initial condition = 1,1
      mutate(invader = "both") %>%
      bind_cols(df.n) %>% 
      bind_rows(df.inv)%>% 
      mutate(external_factor= add_external_factor)
    
    df.sim <- bind_rows(df.sim, df.initc)

    }
  }
}

# save .csv
df.sim <- df.sim %>%
  mutate(function.name = case_when(function.int==1 ~"1.Constant",
                                   function.int==2 ~"2.Linear",
                                   function.int==3 ~"3.Exp",
                                   function.int==4 ~"4.Sigmoid"),
         external_factor = case_when(external_factor=="noise" ~"Noisy change",
                                     external_factor=="none" ~"No external factor",
                                     external_factor=="season" ~"Periodic change"))


write.csv(df.sim , 
          file = paste0("results/df.sim_add_external_factor.csv.gz"))
load("results/df.sim.csv.gz")
df.sim = data.table::fread("results/df.sim_add_external_factor.csv.gz")
head(df.sim)
str(df.sim)
df.sim <- df.sim[,-1]

df.sim$dNi <- as.numeric(df.sim$dNi)
df.sim$dNj <- as.numeric(df.sim$dNj)
df.sim$Ni <- as.numeric(df.sim$Ni)
df.sim$Nj <- as.numeric(df.sim$Nj)


##########################################################################################################
# 2. Visualise population dynamics
##########################################################################################################
sim = 100
t.num= 100
#---- Abundance through time with external factors----
  
example.abundance.external.fact <- df.sim[which(df.sim$sim.i == 2 &
                                                df.sim$time < 51 &
                                                df.sim$invader == "both"),] %>%
  gather(Ni, Nj, key=species, value=abundance) %>%
  ggplot(aes(x=time)) +
  geom_line(aes(y= abundance,color=species)) +
  geom_line(aes(y= external.fact),color="black") +
  theme_bw() + 
  scale_color_manual(values=c("#332288", "#999933")) +
  scale_fill_manual(values=c("#332288", "#999933")) +
  facet_wrap(external_factor ~ function.name ,nrow=3,
             scales="free")+
  labs(title = "Abundance over time of a two-species communities, for one set of parameters,\nfitted in each interaction function respectively",
       y="Local abundance",
       x="Time")+
  theme( legend.key.size = unit(1, 'cm'),
         strip.text = element_text(size=20),
         legend.text=element_text(size=20),
         legend.title=element_text(size=20),
         axis.text.x= element_text(size=16),
         axis.text.y= element_text(size=16),
         #axis.title.y= element_text(size=16),
         title=element_text(size=25))
example.abundance.external.fact
ggsave(example.abundance.external.fact,
       file = "figures/example.abundance.external.fact.pdf")


#---- Abundance through time no external factors , two simulations ----

example.abundance.runawaysim.list <- list()
for(n in c(4,17,67)){
example.abundance.runawaysim.list[[as.character(n)]] <- df.sim[which((df.sim$sim.i == n) &
                                                  df.sim$time < 51 &
                                                  df.sim$external_factor == "No external factor" &
                                                  df.sim$invader == "both"),] %>%
  gather(Ni, Nj, key=species, value=abundance) %>%
  ggplot(aes(x=time)) +
  geom_line(aes(y= abundance,color=species)) +
  geom_line(aes(y= external.fact),color="black") +
  theme_bw() + 
  facet_wrap(.~function.name,nrow=1,scales="free")+
  scale_color_manual(values=c("#332288", "#999933")) +
  scale_fill_manual(values=c("#332288", "#999933")) +
  labs(y="",
       x="") +
  theme( legend.key.size = unit(1, 'cm'),
         strip.text = element_text(size=20),
         legend.text=element_text(size=20),
         legend.title=element_text(size=20),
         axis.text.x= element_text(size=16),
         axis.text.y= element_text(size=16),
         #axis.title.y= element_text(size=16),
         title=element_text(size=25))
}
example.abundance.runawaysim

#######Graph A&B##################

example.abundance.runawaysim <-  ggpubr::ggarrange(plotlist = example.abundance.runawaysim.list, 
                                          ncol=1, nrow=3,
                                          labels=c("Runaway dynamics", 
                                                   "One-species community",
                                                   "Two-species community"),
                            align = "v",
                            hjust=-0.3,
                            vjust = -0.5,    
                            font.label = list(size = 20, color = "black", 
                                              face = "bold", family = NULL, position = "top"),
                            common.legend=T, legend="right") +
  theme(plot.margin = margin(1.5,0.1,0.2,0.1, "cm")) 

library(grid)
example.abundance.runawaysim <- annotate_figure(example.abundance.runawaysim, 
                                                top= textGrob("Abundance over time of a two-species communities, for one set of parameters,\nfitted in each interaction function respectively",  
                                                              gp = gpar(fontsize=20, cex = 1.3)),
                                                left = textGrob("Local abundance", 
                                                                rot = 90, vjust = 3, 
                                                                gp = gpar(fontsize=18,cex = 1.3)),
                                  bottom = textGrob("Time", hjust = 0,vjust=-1.5, 
                                                    gp = gpar(fontsize=18,cex = 1.3)))
example.abundance.runawaysim

ggsave(example.abundance.runawaysim,
       file = "figures/example.abundance.runawaysim.pdf")




#---- Abundance through time for all scenario----
df.min.abundance[which(df.min.abundance$sim.i ==sim &
                         df.min.abundance$function.int == 4  &
                         df.min.abundance$external_factor  =="No external factor"),] 

example.abundance.scenarios <- df.sim[which(df.sim$sim.i ==sim &
                                                       df.sim$function.int == 4  &
                                              df.sim$external_factor  =="No external factor"),] %>%  
  gather(Ni, Nj, key=species, value=abundance) %>%
  ggplot(aes(y=abundance, x=time)) + #geom_point() +
  geom_line(aes(y= abundance,color=species), alpha=0.7) +
  xlim(c(0,100)) + ylim(0,6)+
  geom_smooth(formula = y~x,
              aes(y=abundance, group=species),se=F,
              alpha=0.2,linewidth=0.5, linetype="dashed",color="black") +
  theme_bw() + 
  facet_wrap(as.factor(invader)~.) + 
  scale_color_manual(values=c("#332288", "#999933")) +
  scale_fill_manual(values=c("#332288", "#999933")) +
  labs(title="abundances Ni and Nj over time",
       subtitle=paste0("Growth rate when i is at minimum abundance: ",round(log(GRWL[1,1]),digits = 1),
                      "\nGrowth rate when j is at minimum abundance: ",round(log(GRWL[1,2]),digits = 1))) 
example.abundance.scenarios
ggsave(example.abundance.scenarios,
       file = "figures/example.abundance.scenarios.pdf")


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



