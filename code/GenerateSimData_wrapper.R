# Generate Data based on Stouffer 2022 - Methods in Ecology and Evolution

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation of 100 pairwise communities with the same parameters for one season
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# the community dynamics over 10 seasons is as followed:
source("code\GenerateConitunousData-Stouffer.R")

# communities fecundities and abundances are stored on two csv
Generate.simulated.data <- NULL
Generate.experimental.outcomes <- NULL
for(ri in c(1:100)){
  print(ri)
random.set <- ri
set.seed(random.set) #to create reproducible results
source("code/GenerateSimData-Stouffer.R") # Generate Simulated data in "simulated.data" df and 
simulated.data$sim <- ri
Generate.simulated.data <- bind_rows(Generate.simulated.data,
                                     simulated.data)
experimental.outcomes$sim <- ri
Generate.experimental.outcomes <- bind_rows(Generate.experimental.outcomes,
                                            experimental.outcomes)
}

write_csv(Generate.experimental.outcomes, 
          file = "results/Generate.experimental.outcomes.csv")
write_csv(Generate.simulated.data, 
          file = "results/Generate.simulated.data.csv")

# visualisation of communities fecundity distributions


ggsave("figures/simulated.seed.density.pdf",
       plot = ggplot(Generate.simulated.data,
                     aes(x=fecundity,
                         group = sim)) + 
         geom_density(alpha=0.6) +
         #scale_color_manual(values=c("blue","red")) +
         facet_wrap(focal~.) +
         xlab("Viable seed, fecundity per individuals") + 
         theme_bw()
)
ggsave("figures/simulated.seed.biomass.pdf",
       plot = ggplot(Generate.experimental.outcomes) + 
         geom_point( aes(x=biomass.i,y=biomass.j,
                         color = sim),alpha=0.3) +
         #geom_smooth( aes(x=biomass.i,y=biomass.j),color="grey",alpha=0.6) +
         #xlim(0.875,1)+
         #ylim(0,0.6) +
         theme_bw()+
         labs(title="Biomass of species j in function of biomass of species i")
)


#---- others -----
#### MALYON's model ####
set.seed(1236) #to create reproducible results
source("code/GenerateSimData-Malyon.R")
simdata$fecundity <- as.numeric(simdata$fecundity)
simdata$Spi <- as.numeric(simdata$Spi)
simdata$Spj <- as.numeric(simdata$Spj)
ggsave("figures/simulated.seed.density.pdf",
       plot = ggplot(simdata, aes(x=seeds, fill=as.factor(focal), color=as.factor(focal))) + 
         geom_density(alpha=0.5) + 
         scale_x_continuous(trans='log2',
                            name = "Viable seed, fecundity per individuals (ln transformed)") + 
         theme_bw() + scale_color_manual(values=c("blue","red")) + 
         scale_fill_manual(values=c("blue","red"))
)

#### Ricker model Data ####
source("code/simul_data.R")
set.seed(16)
simul_data_ricker <- simul_data(S = 2,   # number of focal groups / species
                                K = 2,   # number of neighbour focals 
                                pI = 0.1) # Generate Simulated data in "simulated.data" df and
simdata <- simul_data_ricker$simdata
cols.num <- c("seeds","i","j")
simdata[cols.num] <- sapply(simdata[cols.num],as.numeric)

ggsave("figures/Ricker.fecundity.pdf",
       plot = ggplot(simdata, aes(x=log(seeds), color=as.factor(focal))) + 
         geom_density(alpha=0.5) + 
         scale_x_continuous(name = "Viable seed, fecundity per individuals (log transformed)") + 
         theme_bw() 
)
simul_data_ricker$sim_alpha_specific


