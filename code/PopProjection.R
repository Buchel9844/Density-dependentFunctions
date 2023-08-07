# Project population trajectory of species i and j for each scenario and each function ? 

library(deSolve)

#----extract parameters for each function and each scenario----
params.low <- list(
  g = c(0.7, 0.7), # germination rate of seeds
  s    = c(0.9, 0.95) # 1 - mortality rate of seeds
)
params.medium <- list(
  g = c(0.7, 0.7), # germination rate of seeds
  s    = c(0.90, 9) # mortality rate of seeds
)
params.high <- list(
  g = c(0.7, 0.7), # germination rate of seeds
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
      simdata <- simdata[which(simdata$time >= time.exp.min & 
                                 simdata$time < time.exp.max &
                                 simdata$focal == Code.focal),]
      Nmax[int.focal,1] <- simdata[which.max(simdata$fecundity),"plants.i"]
      Nmax[int.focal,2] <- simdata[which.max(simdata$fecundity),"plants.j"]
      
      
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
population.dynamics.all <- data.frame()
for( scenario in c("low","medium","high")){
population.dynamics <- read.csv(paste0("results/Generate.population.dynamics.",scenario,".csv"))
population.dynamics$scenario <- scenario

# initial species densities
state <- c(population.dynamics[which(population.dynamics$Time == time.exp.max),"Seeds.i"][1],
           population.dynamics[which(population.dynamics$Time == time.exp.max),"Seeds.j"][1]) 

state.list[[scenario]] <- state
population.dynamics.all  <- bind_rows(population.dynamics.all ,population.dynamics)
}
population.dynamics.all <- pivot_longer(population.dynamics.all , 
                              cols=c("Seeds.i","Seeds.j"), 
                              names_to="Population", values_to="N") %>%
  mutate(Population = case_when(Population =="Seeds.i" ~"Ni" ,
                                Population =="Seeds.j" ~"Nj" ))
write.csv(population.dynamics.all, 
          "results/population.dynamics.all.csv")


#---- for high scenario ----
source("code/PopProjection_toolbox.R")
df.pop.proj <- NULL

for( scenario in c("low","medium","high")){
  for(function.int in c(1:4)){ # c(1:4)
    df <- NULL
p <- list.pars[[paste0("pars_",scenario,"_function_",function.int)]] 
y <- state.list[[scenario]]  # initial species densities
gens = 25

df <- Ricker_solution(gens= gens,
                state=y,
                pars=p)
df$scenario <- scenario
df$function.int <- function.int
df.pop.proj <- bind_rows(df.pop.proj,df)
  }
}



df.pop.proj.l <- pivot_longer(df.pop.proj , cols=2:3, names_to="Population", values_to="N")
df.pop.proj.l <- df.pop.proj.l %>% arrange(factor(scenario, 
                                                  levels = c("low","medium","high")))
write.csv(df.pop.proj.l, 
          "results/PopProjection.csv")

ggplot() + geom_line(data=df.pop.proj.l,aes(x=t, y=log10(N), colour=Population)) +
  geom_point(data=df.pop.proj.l,aes(x=t, y=log10(N), colour=Population)) +
  facet_grid(scenario~function.int, scale="free")+
  geom_line(data=population.dynamics.all,
            aes(x=Time-time.exp.max, y=log10(N), colour=Population),
            alpha=0.5) +
  theme_bw() +
  scale_color_manual(values=c("black","red")) +
  
  labs(title = "Projection of seed trajectory for species i (native) and j (competitor)",
       x="Generations", y="Density of viable seeds")


#########################
#Coexistence probability 
########################
df.prob.coexist <- NULL
for( scenario in c("low","medium","high")){
  for(function.int in c(1:4)){ # c(1:4)
    df <- NULL
    p <- list.pars[[paste0("pars_",scenario,"_function_",function.int)]] 
    y <- state.list[[scenario]]  # initial species densities
    gens = 25
    
    df <- grwr(par.dat = p, t.num  = 100)
    df$scenario <- scenario
    df$function.int <- function.int
    
    df.prob.coexist <- bind_rows( df.prob.coexist,df)
  }
}

df.prob.coexist.out <- df.prob.coexist %>%
  dplyr::select(invader, grwr, scenario ,function.int , Cgrwc) %>%
  unique()  %>%
  mutate(grwrChesson = log(grwr)-log(Cgrwc))
  

ggplot(df.prob.coexist.out,aes(x=invader,y=grwrChesson)) +
  geom_point(aes(color=as.factor(function.int)),size=2,alpha=0.6)+
  facet_grid(.~scenario )
