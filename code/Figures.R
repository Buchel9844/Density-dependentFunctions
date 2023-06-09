# This script will make the figures of chap 2, number- XXX of the manuscript

library(ggplot2)
library(ggridges)
library(tidyverse)
##########################################################################################################
# 1. Simulation data
##########################################################################################################
#---- 1.0. Compare simulation with estimate alpha for ricker model  ----

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
load("results/Alphadistribution.neighbours.csv.gz")


AlphadistributionGraph <- ggplot(Alphadistribution.neighbours, aes(x= abundance.neighbours,
                                                           y=alpha_mean)) + 
  geom_smooth(alpha=0.5,se=F, aes(
    color=density.function,
    fill=density.function)) +
  geom_point(aes(
    color=density.function,
    fill=density.function)) +
  #ylim(-0.1,0.1)+
  geom_hline(yintercept=0,color="dark grey") +
  xlim(0,max(Alphadistribution.neighbours$abundance.neighbours)) + 
  geom_errorbar(aes(ymin=alpha_mean-alpha_sd^2, ymax=alpha_mean+alpha_sd^2,
                    color=density.function), width=.2,
                              position=position_dodge(0.05)) +
  facet_grid( focal ~ neighbours,scale="free", switch="both") + theme_bw() + 
  guides(fill="none",color="none")+
  scale_color_manual("Density-dependent functions",values=c("black","#CC79A7","#E69F00","#009E73")) +
  labs(title="Direct interactions distributions of the 4 density-dependent functions", 
       y="Resulting effect", x=" Neighbour density ")

ggsave("figures/AlphadistributionGraph.pdf",
       plot = AlphadistributionGraph
)


#---- 1.2. Figure of fecundity estimation distribution ----

load("results/Fecunditydistribution.csv.gz")
#APM_Fecunditydistribution <- read.csv("results/APM_PostFecunditydistribution.csv")
load("results/PostFecunditydistribution.csv.gz")


#Post_Fecunditydistribution <- bind_rows(Post_Fecunditydistribution,
#                                        APM_Fecunditydistribution)

simdata <- read.csv("results/simulated.data.csv")

Fecunditydistribution.n <- data.frame(seed=simdata$fecundity,
                                    focal=simdata$focal,
                                    obervation= NA )
                                    #max.seed= NA,
                                    #alpha.function = "initial.seed"
                                    # )%>%
 # group_by(focal) %>%
  #mutate(max.seed= max(seed,na.rm=T))

Post_Fecunditydistribution <- PostFecunditydistribution %>%
  filter(!is.na(Fec))

PostFecundityGraph <- ggplot(PostFecunditydistribution) +
  geom_density(aes(x=Fec,
                  group=iterations,
                  color=as.factor(alpha.function)),
               alpha=0.2) +
  geom_density(data = Fecunditydistribution.n , aes(x=seed, 
                                          color=as.factor(focal)),
               alpha=0.5) +
  facet_wrap(focal ~  alpha.function, nrow=2, scales = "free" ) +
  scale_x_continuous(limits = function(x) {x + c(-10, 10)}) +
  labs(title = " Fecundity distributions for predictions compared to initial data", y ="Density",
       x="Fecundity") +  theme_bw() + guides(color="none") + 
  scale_color_manual(values=c("black","#CC79A7","#E69F00","#009E73",
                              "blue","red"))

ggsave("figures/PostFecundityGraph.pdf",
       plot = PostFecundityGraph)


##########################################################################################################
# 2. Empirical data
##########################################################################################################

#---- 2.1. Figure of empirical seed distribution ----

