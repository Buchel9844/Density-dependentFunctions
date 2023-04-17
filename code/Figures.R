# This script will make the figures of chap 2, number- XXX of the manuscript

library(ggplot2)
library(ggridges)
library(tidyverse)
##########################################################################################################
# 1. Simulation data
##########################################################################################################
#---- 1.0. Compare simulation with estimate alpha  ----

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
#---- 1.1. Figure of density-dependent function ----
Alphadistribution.neighbours <- read.csv("results/Alphadistribution.neighbours.csv")


AlphadistributionGraph <- ggplot(Alphadistribution.neighbours, aes(x= abundance.neighbours,
                                                           y=alpha_mean)) + 
  geom_smooth(alpha=0.5,se=F, aes(
    color=density.function,
    fill=density.function)) +
  geom_point(aes(
    color=density.function,
    fill=density.function)) +
  ylim(-5,5)+
  xlim(0,max(Alphadistribution.neighbours$abundance.neighbours)) + 
  geom_errorbar(aes(ymin=alpha_mean-alpha_sd, ymax=alpha_mean+alpha_sd,
                    color=density.function), width=.2,
                              position=position_dodge(0.05)) +
  facet_grid(neighbours ~ focal,scale="free") + theme_bw()

ggsave("figures/AlphadistributionGraph.pdf",
       plot = AlphadistributionGraph
)


#---- 1.2. Figure of fecundity estimation distribution ----

Fecunditydistribution <- read.csv("results/Fecunditydistribution.csv")
#APM_Fecunditydistribution <- read.csv("results/APM_PostFecunditydistribution.csv")
Post_Fecunditydistribution <- read.csv("results/PostFecunditydistribution .csv")


#Post_Fecunditydistribution <- bind_rows(Post_Fecunditydistribution,
#                                        APM_Fecunditydistribution)

#simulated.data <- read.csv("results/simulated.data.csv")

Fecunditydistribution.n <-data.frame(seed=simdata$seeds,
                                    focal=paste("species",simdata$focal),
                                    obervation= NA, 
                                    max.seed= NA,
                                    alpha.function = "initial.seed"
                                     )%>%
  group_by(focal) %>%
  mutate(max.seed= max(seed,na.rm=T))

Post_Fecunditydistribution <- Post_Fecunditydistribution%>%
  filter(!is.na(Fec))%>%
  group_by(focal) %>%
  mutate(max.seed = max(simdata$seeds,na.rm=T)+20)%>%
  filter(Fec < max.seed)


PostFecundityGraph <- ggplot(Post_Fecunditydistribution) +
  geom_density(aes(x=Fec,y=after_stat(scaled),
                                               group=iterations),
               alpha=0.5,color="grey") + theme_bw() +
  facet_wrap(focal ~  alpha.function, nrow=2, scales = "free") +
  geom_density(data = simdata, aes(x=seeds, 
                                          y=after_stat(scaled)),
               alpha=0.5,color="black") 

ggsave("figures/PostFecundityGraph.pdf",
       plot = PostFecundityGraph)

ggplot(Fecunditydistribution.n,aes(seed,y=after_stat(scaled))) + geom_density(color="grey")+
  geom_density(data = simulated.data, aes(x=fecundity, 
                                          y=after_stat(scaled)),
               alpha=0.5,color="black") 

##########################################################################################################
# 2. Empirical data
##########################################################################################################

#---- 2.1. Figure of empirical seed distribution ----

