# Generate Data based on Stouffer 2022 - Methods in Ecology and Evolution

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation of 100 pairwise communities with the same parameters for one season
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(Stouffer == T){
# communities fecundities and abundances are stored on two csv
  
random.set <- 1616
set.seed(random.set) #to create reproducible results
source("code/GenerateSimData-Stouffer.R") # Generate Simulated data in "simulated.data" df and 

write_csv(experimental.outcomes, 
          file = paste0("results/Generate.experimental.outcomes.",scenario,".csv"))
write_csv(simulated.data, 
          file = paste0("results/Generate.simulated.data.",scenario,".csv"))

# visualisation of communities fecundity distributions

# the community dynamics over 10 seasons is as followed:
source("code/GenerateContinuousData-Stouffer.R")


ggsave(paste0("figures/simulated.seed.density.",scenario,".pdf"),
       plot = ggplot(simulated.data[which(simulated.data$time>=time.exp.min & simulated.data$time < time.exp.max),],
                     aes(x=fecundity,
                         group=focal)) + 
         geom_density(alpha=0.6) +
         #scale_color_manual(values=c("blue","red")) +
         labs(x=paste0("Viable seed, fecundity per individuals for scenario ",scenario),
             title=paste0("Seed distribution (exp), for both individuals for scenario ",scenario)) + 
         theme_bw()
)

ggsave(paste0("figures/simulated.seed.biomass.",scenario,".pdf"),
       plot = ggplot(experimental.outcomes) + 
         geom_point( aes(x=biomass.i,y=biomass.j),alpha=0.3) +
         #geom_smooth( aes(x=biomass.i,y=biomass.j),color="grey",alpha=0.6) +
         #xlim(0.875,1)+
         #ylim(0,0.6) +
         theme_bw()+
         labs(title=paste0("Biomass of species j in function of biomass of species i for scenario ",scenario))
)


ggsave(paste0("figures/simulated.seed.biomass.",scenario,".pdf"),
       plot = ggplot(data = gather(simulated.data[which( simulated.data$time>=time.exp.min & simulated.data$time < time.exp.max),],
                           plants.i,plants.j,key="plants",value="n.ind"),
                    aes(y=fecundity, x=time, color=n.ind)) +
                      geom_point(aes(shape=focal),size=3,position="jitter")+ theme_bw() +
           labs(title=paste0("Fecundity of species j in the first ",time.exp," for scenario ", scenario))
)

}
#---- others -----
#### MALYON's model ####
if(Malyon == T){
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
}
#### Ricker model Data ####
if(Ricker == T){
source("code/GenerateSimData_Ricker.R")
set.seed(16)
simul_data_ricker <- simul_data(S = 2,   # number of focal groups / species
                                K = 2,   # number of neighbour focals 
                                pI = 0.1) # Generate Simulated data in "simulated.data" df and
simdata <- simul_data_ricker$simdata
cols.num <- c("seeds","i","j")
simdata[cols.num] <- sapply(simdata[cols.num],as.numeric)

ggsave("figures/Ricker.fecundity.pdf",
       plot = ggplot(simdata, aes(x=seeds, color=as.factor(focal))) + 
         geom_density(alpha=0.5) + 
         scale_x_continuous(name = "Viable seed, fecundity per individuals (log transformed)") + 
         theme_bw() 
)
simul_data_ricker$sim_alpha_specific
}

