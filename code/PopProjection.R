# Project population trajectory of species i and j for each scenario and each function ? 

library(deSolve)

#----extract parameters for each function and each scenario----
params.low <- list(
  g = c(0.001, 0.001), # germination rate of seeds
  s    = c(0.9, 0.95) # 1 - mortality rate of seeds
)

params.medium <- list(
  g = c(0.001, 0.001), # germination rate of seeds
  s    = c(0.90, 0.85) # mortality rate of seeds
)
params.high <- list(
  g = c(0.001, 0.001), # germination rate of seeds
  s    = c(0.90, 0.50) # mortality rate of seeds
)

list.pars <- list()
for( scenario in c("low","medium","high")){

   for (function.int in c(1:4)){ # c(1:4)
     a_initial <-matrix(nrow=2,ncol=2)
     a_slope <-matrix(nrow=2,ncol=2)
     c <- matrix(nrow=2,ncol=2) 
     Nmax <- matrix(nrow=2,ncol=2)  
     lambda <-c()
      for(Code.focal in c("i","j")){ #,"j"
      print(paste(scenario, Code.focal,", function",function.int))
      
      load(paste0("results/FinalFit_",scenario,"_",Code.focal,"_function_",function.int,".rds"))

      if(Code.focal =="i"){int.focal <- 1}else{int.focal <- 2}
      lambda[int.focal] <- c(summary(FinalFit)$summary["lambdas[1]","mean"])
      
      a_initial[int.focal,] <- c(summary(FinalFit)$summary["alpha_initial[1]","mean"],
                                 summary(FinalFit)$summary["alpha_initial[2]","mean"])
      a_slope[int.focal,] <- c(summary(FinalFit)$summary["alpha_slope[1]","mean"],
                                 summary(FinalFit)$summary["alpha_slope[2]","mean"])
      c[int.focal,] <- c(summary(FinalFit)$summary["c[1]","mean"],
                               summary(FinalFit)$summary["c[2]","mean"])
  
      
      simdata <- read.csv(paste0("results/Generate.simulated.data.",scenario,".csv"))
      simdata <- simdata[which(simdata$time <= time.exp &
                                 simdata$focal == Code.focal),]
      Nmax[int.focal,1] <- 0 #simdata[which.max(simdata$fecundity),"plants.i"]
      Nmax[int.focal,2] <- 0 #simdata[which.max(simdata$fecundity),"plants.j"]
      
      
      remove(FinalFit )

      }
     #population.dynamics <- read.csv(paste0("results/Generate.population.dynamics.",scenario,".csv"))
     # initial species densities
    #Nmax[1] <- population.dynamics[which(population.dynamics$Time == 2),"Seeds.i"][1]
     #Nmax[2] <- population.dynamics[which(population.dynamics$Time == 2),"Seeds.j"][1]
     
     list.pars[[paste0("pars_",scenario,"_function_",function.int)]] <- append( list(lambda = lambda, 
                                                                                     a_slope = a_slope,
                                                                                     a_initial =  a_initial,
                                                                                     c=c,
                                                                                     Nmax = Nmax,
                                                                                     function.int = function.int ),
                                                                                     get(paste0("params.",scenario)))
       
   }
}
state.list <- list()
for( scenario in c("low","medium","high")){
population.dynamics <- read.csv(paste0("results/Generate.population.dynamics.",scenario,".csv"))
# initial species densities
state <- c(population.dynamics[which(population.dynamics$Time == 3),"Seeds.i"][1],
           population.dynamics[which(population.dynamics$Time == 3),"Seeds.j"][1]) 

state.list[[scenario]] <- state
}

#---- for high scenario ----
source("code/PopProjection_toolbox.R")
df.pop.proj <- NULL

for( scenario in c("low","medium","high")){
  for(function.int in c(1:4)){ # c(1:4)
    df <- NULL
p <- list.pars[[paste0("pars_",scenario,"_function_",function.int)]] 
y <- state.list[[scenario]]  # initial species densities
gens= 20

df <- Ricker_solution(gens= 20,
                state=state,
                pars=pars)
df$scenario <- scenario
df$function.int <- function.int
df.pop.proj <- bind_rows(df.pop.proj,df)
  }
}

df.pop.proj.l <- pivot_longer(df.pop.proj , cols=2:3, names_to="Population", values_to="N")
ggplot(aes(t, N, colour=Population), data=df.pop.proj.l) + geom_line() +
  facet_grid(function.int ~scenario, scale="free")


years = seq(0, 20, by = 1)

df.NiNj <- ode(y = y, times = years, 
                             func = Ricker_function, parms = p,
               method = "euler")


plot(df.NiNj)
head(df.NiNj)
df.NiNj.l <- pivot_longer(as.data.frame(df.NiNj), 
                          cols=2:3, names_to="Population", values_to="N")

ggplot(aes(time, N, colour=Population), data=df.NiNj.l) + geom_line() +
  scale_y_log10()


df.NiNj$function.int <- function.int
df.NiNj$scenario <- scenario

df.pop.proj <- bind_rows(df.pop.proj, df.NiNj)

