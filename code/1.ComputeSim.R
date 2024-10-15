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
library(ggplotify) 
library(odeintr)
library(brms)
library(bbmle)
library(reshape2)
library(ggthemes)
library(tsvr) # library for the timescale specific variance ratio
library(data.table) # write csv.gz
library(ppcor) # partial correlation value
library(scales) # for pretty breaks with scale_axis
library(grid)

###########################################################################################################
# 1. # 1. Compute Low- Density Growth Rates for wide range of parameters without coexistence restriction
##########################################################################################################
source("code/1.1.PopProjection_toolbox.R")
# define parameters
nsims <- 1500
species <- 2
function.int <- c(1:4)
t.num = 200
gen = seq(from=0, to=t.num , by=1)
# run sim


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

white.noise <- function(t.num, b){
  e_noise = rnorm(t.num+1,mean=0, sd= b) 
  for( n in 1:sample(2:10)[1]){
    e_noise[n] <- rnorm(1,mean=0,
                        sd=abs(rnorm(1, mean=0, sd=b*5)))
  }
  return( e_noise)
}

set.seed(1644)

params <-list()
i <- 1
init.facilitation_tally <- 0
init.comp_tally <- 0
init.faccomp_tally <- 0
for( n in c(1:20000)){
  if(i==nsims+1) next
  print(paste(n,i,init.facilitation_tally,init.comp_tally,init.faccomp_tally))
  b <-  sample(seq(0,2,0.01),1) # amplitude of external factor
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
                    
                    a_initial = matrix(ncol=species, nrow = species, 
                                       data = rtruncnorm(n=species*species, mean= -0.1, sd=0.5, a=-1, b=1),
                                       dimnames = list(c("i", "j"),
                                                       c("i", "j"))),
              
                    c =matrix(ncol=species, nrow = species, 
                                     data =runif(n=species*species, min=-1, max=0),
                              dimnames = list(c("i", "j"),
                                              c("i", "j"))),
                    b =  b,
                    e_seasonal = sin((2*pi/20)*gen)*abs(rnorm(1,mean=0, sd=b)),
                    e_noise = white.noise(t.num, b)
                  
                    )
    
    params[[i]]$a_initial[1,1] <- -abs(params[[i]]$a_initial[1,1])
    params[[i]]$a_initial[2,2] <- -abs(params[[i]]$a_initial[2,2])
    
    
    if(params[[i]]$a_initial[1,2] > 0 & params[[i]]$a_initial[2,1] > 0){
      if(init.facilitation_tally >= nsims/3) {
        i = i 
        next
      }else{i = i + 1
      init.facilitation_tally <- init.facilitation_tally + 1
      next
      } 
    }
    if((params[[i]]$a_initial[1,2] > 0 & params[[i]]$a_initial[2,1] < 0 )|
       (params[[i]]$a_initial[1,2] < 0 & params[[i]]$a_initial[2,1] > 0)){
      if(init.faccomp_tally >= nsims/3){
        i = i
        next
      }else{
        i = i + 1
        init.faccomp_tally <- init.faccomp_tally + 1
        next
      }
    }
        
    if(params[[i]]$a_initial[1,2] < 0 & params[[i]]$a_initial[2,1] < 0){
      if(init.comp_tally >= nsims/3){
        i = i
        next
      }else{
        i = i + 1
        init.comp_tally <- init.comp_tally + 1
        next
         }
       }
      
}

#check all tally equal nsims/3
length(params)
init.facilitation_tally ==nsims/3
init.comp_tally==nsims/3
init.faccomp_tally ==nsims/3



save(params,
     file = paste0("results/sim_params.Rdata"))
load("results/sim_params.Rdata")

df.sim  <- NULL

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
    
    #df.inv <- NULL
    #df.inv <- GrowthSimInv(par.dat = p, t.num  = t.num,
    #     function.int,
    #    add_external_factor) %>%# simulation of invasion for both species
    #  bind_cols(df.n) 
   
    df.initc <- Ricker_solution_withalpha(state = c(5,5), pars= p, gens=t.num,
                        function.int,add_external_factor) %>% # simulation of initial condition = 1,1
      mutate(invader = "both") %>%
      bind_cols(df.n) %>% 
     # bind_rows(df.inv)%>% 
      mutate(external_factor= add_external_factor)
    
    df.sim <- bind_rows(df.sim, df.initc)

    }
  }
}

# save .csv
df.sim <- df.sim %>%
  mutate(function.name = case_when(function.int==1 ~"1.Traditional",
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
  
example.abundance.external.fact <- df.sim[which(df.sim$sim.i ==500 &
                                                df.sim$time < 100 &
                                                df.sim$invader == "both"),] %>%
  gather(Ni, Nj, key=species, value=abundance) %>%
  mutate(abundance = case_when(abundance > 5000 ~ 5000, 
                               T ~ abundance)) %>%
  ggplot(aes(x=time)) +
  geom_line(aes(y= log(abundance+1),color=species)) +
  #geom_line(aes(y= external.fact),color="black") +
  theme_bw() + 
  scale_color_manual(values=c("#332288", "#999933")) +
  scale_fill_manual(values=c("#332288", "#999933")) +
  facet_wrap(external_factor ~ function.name ,nrow=3)+
  labs(title = "Abundance over time of a two-species communities, for one set of parameters,\nfitted in each functional form respectively",
       y="Local abundance (log)",
       x="Time")+
  scale_y_continuous(expand=c(0.01,0.01)) +
  theme( legend.key.size = unit(2, 'cm'),
         strip.text = element_text(size=20),
         strip.background = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor = element_blank(),
         legend.text=element_text(size=20),
         legend.title=element_text(size=20),
         axis.text.x= element_text(size=16),
         axis.text.y= element_text(size=16),
         #axis.title.y= element_text(size=16),
         title=element_text(size=25))
example.abundance.external.fact
 ggsave(example.abundance.external.fact,
       file = "figures/example.abundance.external.fact.specific.case.pdf")


#---- Abundance through time no external factors , two simulations ----

# Facilitation 
comp <- df.sim$sim.i[which(df.sim$a_initial.i.j < -0.1 &
                                df.sim$a_initial.i.j > df.sim$a_initial.i.i &
                                df.sim$a_initial.j.i > df.sim$a_initial.j.j &
               df.sim$a_slope.i.j < -0.1 &
               df.sim$c.j.i< -0.7 &
               df.sim$time == 1 &
               df.sim$function.int == 1 &
               df.sim$external_factor == "No external factor" &
               df.sim$invader == "both")]

 fac <- df.sim$sim.i[which(df.sim$a_initial.i.j > 0 &
                             df.sim$time == 1 &
                             df.sim$function.int == 1 &
                             df.sim$function.int == 1 &
                             df.sim$external_factor == "No external factor" &
                             df.sim$invader == "both")]
 
 fac <- df.sim$sim.i[which(df.sim$a_initial.i.j > 0 &
                             df.sim$a_slope.i.j < - 0.15 &
                             df.sim$a_slope.i.j > - 0.3 &
                             df.sim$c.i.j< -0.7 &
                             df.sim$time == 1 &
                             df.sim$function.int == 1 &
                             df.sim$external_factor == "No external factor" &
                             df.sim$invader == "both")]

roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

fancy_scientific <- function(l){
  l <- roundUpNice(l)
  #l <- round(l)
  # turn in to character string in scientific notation
  l <- format(l)
  # round to one point after the .
  l <- gsub("^(.*)e", "'\\1'e", l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  return(l)
}
comp
fac  
my_cols <- c("#AA4499", "#DDCC77","#88CCEE", "#44AA99")

#584,156
#280
params[[490]]$a_initial
{
example.abundance.runawaysim.list <- list()
 for(n in  c(490,28)){
   df.plot <- NULL
   df.alpha <-NULL
    for( fct.int in 1:4){
      df.plot.n <- Ricker_solution_withalpha(gens=21,
                                             state=c(5,5),pars = params[[n]],
                                             function.int=fct.int,
                            add_external_factor="none") %>%
        add_column("function.int"= fct.int) %>%
        filter(time < 21)

      df.plot.n$resulting.alpha.i <-  df.plot.n$aij*df.plot.n$Nj + df.plot.n$aii*df.plot.n$Ni 
      df.plot.n$resulting.alpha.j <-  df.plot.n$aji*df.plot.n$Ni + df.plot.n$ajj*df.plot.n$Nj 
      df.plot <- bind_rows(df.plot, df.plot.n)
      
      
      df.alpha_n <- as.data.frame(params[[n]][c("sim","Nmax","a_slope","a_initial","c")]) %>%
        mutate(focal = c("i","j")) %>%
        mutate(function.int = fct.int ) %>%
        filter(focal =="i") %>%
        remove_rownames() %>%
        bind_cols(data.frame(density=c(0:20)))%>%
        as.data.frame()
      
        if(df.alpha_n$function.int[1]==1){
          df.alpha_n$alpha_inter_value <- df.alpha_n$a_initial.j[1]
          
        }
        if(df.alpha_n$function.int[1]==2){
          df.alpha_n$alpha_inter_value <- alpha_function2(df.alpha_n$a_initial.j[1],
                                                          df.alpha_n$a_slope.j[1],
                                                          df.alpha_n$density,
                                                          df.alpha_n$Nmax.j[1])
          
        }
        if(df.alpha_n$function.int[1]==3){
          df.alpha_n$alpha_inter_value <- alpha_function3(df.alpha_n$a_initial.j[1],
                                                          df.alpha_n$a_slope.j[1],
                                                          df.alpha_n$c.j[1],
                                                          df.alpha_n$density,
                                                          df.alpha_n$Nmax.j[1])
          
        }
        if(df.alpha_n$function.int[1]==4){
          df.alpha_n$alpha_inter_value <- alpha_function4(df.alpha_n$a_initial.j[1],
                                                          df.alpha_n$a_slope.j[1],
                                                          df.alpha_n$c.j[1],
                                                          df.alpha_n$density,
                                                          df.alpha_n$Nmax.j[1])
        }
        
        df.alpha <- bind_rows(df.alpha,df.alpha_n )
        
      }

   df.plot <- df.plot %>% mutate(function.name = case_when(function.int==1 ~"Constant",
                                    function.int==2 ~"2.Linear",
                                    function.int==3 ~"3.Exp",
                                    function.int==4 ~"4.Sigmoid"))
   
   df.alpha <-  df.alpha %>% mutate(function.name = case_when(function.int==1 ~"1.Constant",
                                                           function.int==2 ~"2.Linear",
                                                           function.int==3 ~"3.Exp",
                                                           function.int==4 ~"4.Sigmoid"))
   
#scale.resulting.alpha <- c(df.plot[,c("resulting.alpha.i")],df.plot[,c("resulting.alpha.j")])
scale.resulting.alpha <- c(df.plot[,c("resulting.alpha.i")])

scale.resulting.alpha[which(scale.resulting.alpha < -1)] <- -1
scale.resulting.alpha[which(scale.resulting.alpha > 1)] <- 1

plot.abundance <- df.plot %>%
  gather(Ni, Nj, key=species, value=abundance) %>%
  filter(species=="Ni") %>%
  #mutate(abundance = case_when(is.na(abundance) ~ 0 ,
  #                             abundance > 10000 ~ 10000,
   #                            T ~ abundance)) %>%
  add_column(#"abundance.label" = abundance,
    scale.resulting.alpha= scale.resulting.alpha) %>%
  ggplot(aes(x=time)) +
  geom_line(aes(y= log(abundance+1),color=function.name),size=3,alpha=0.8) +
  theme_bw() + 
  labs(x="Time in years",y="Local abundance of\nfocal species (log)", #y="",#
       color= "") +
  scale_color_manual(values=my_cols, labels=c("Traditional\nconstant",
                                              "Linear",
                                              "Exponential",
                                              "Sigmoid")) +
  coord_cartesian(ylim=c(0,8)) +
  #guides(color="none") +
  theme( legend.key.size = unit(2, 'cm'),
         strip.text = element_text(size=20),
         strip.background = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor = element_blank(),
         legend.text=element_text(size=20),
         legend.title=element_text(size=20),
         axis.text.x= element_text(size=14),
         axis.text.y= element_text(size=14),
         #axis.title.y= element_text(size=16),
         title=element_text(size=25))

plot.alpha <- df.alpha %>%
  ggplot(aes(x=density)) +
  geom_point(aes(y= alpha_inter_value,color=function.name),size=3,alpha=0.8) +
  geom_line(aes(y= alpha_inter_value,color=function.name),size=3,alpha=0.8) +
  theme_bw() + 
  geom_hline(yintercept = 0) + 
  scale_y_continuous( expand=c(0,0)) +
  coord_cartesian(ylim=c(-1.5,1.5)) +
  geom_label(x=16,y=1.35,label="Facilitation",size=8) +
  geom_label(x=16,y=-1.35,label="Competition",size=8) +
  scale_color_manual(values=my_cols,labels=c("Traditional\nconstant",
                                            "Linear",
                                            "Exponential",
                                            "Sigmoid")) +
  labs(color="", y="Interspecific\nPer capita effect" ,#y="",#
       #title=n,
       x="Neigbours density") +
  theme( legend.key.size = unit(2, 'cm'),
         legend.position = "bottom",
         strip.text = element_text(size=20),
         strip.background = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor = element_blank(),
         legend.text=element_text(size=20),
         legend.title=element_text(size=20),
         axis.text.x= element_text(size=14),
         axis.text.y= element_text(size=14),
         #axis.title.y= element_text(size=16),
         title=element_text(size=25))


plot.legend <- get_legend(plot.alpha)
  
if(length(example.abundance.runawaysim.list) ==0){
  example.abundance.runawaysim.list[[paste0("plot",n)]]  <- ggarrange(plot.alpha,ggplot() + theme_void(),plot.abundance,
                                                                      ncol=3,nrow=1,
                                                                      widths = c(1,0.1,1),
                                                                      common.legend = T,
                                                                      legend="none")
}
if(length(example.abundance.runawaysim.list) ==1){
  example.abundance.runawaysim.list[[paste0("plot",2)]]  <- ggplot() + theme_void()
  
}
if(length(example.abundance.runawaysim.list) >= 2){
  example.abundance.runawaysim.list[[paste0("plot",n)]]  <- ggarrange(plot.alpha,ggplot() + theme_void(),plot.abundance,
                                                     ncol=3,nrow=1,
                                                     widths =c(1,0.1,1),
                                                     common.legend = T,
                                                     legend="none")

 }
 }

  
example.abundance.runawaysim.list[["legend"]] <- ggplotify::as.ggplot(plot.legend)
example.abundance.runawaysim.list$plot490
length(example.abundance.runawaysim.list)
#######Graph A&B##################

example.abundance.runawaysim <-  ggpubr::ggarrange(plotlist = example.abundance.runawaysim.list, 
                                           nrow=4, ncol=1,
                                           heights = c(1,0.2,1,0.2),
                                          labels=c("a. Predominant interaction is facilitative", 
                                                   "",
                                                   "b. Predominant interaction is competitive",
                                                   ""),
                            align = "v",
                            hjust= 0,
                            vjust = -0.5,    
                            font.label = list(size = 28, color = "black", 
                                              face = "bold", family = NULL, position = "top"),
                           common.legend = T, legend="bottom") +
  theme(plot.margin = margin(1.3,1,0,1, "cm"))  
 example.abundance.runawaysim
}
library(grid)
 
example.abundance.runawaysim <- annotate_figure(#example.abundance.runawaysim, left = textGrob("Local abundance (log)", 
                                                  #              rot = 90, vjust = 3, 
                                                  #              gp = gpar(fontsize=18,cex = 1.3)),
                                              #  right = textGrob("Facilitation", 
                                                                # check.overlap = T,
                                                  #               rot=90,
                                                   #             vjust = -8, hjust = 0.25,
                                                     #           gp = gpar(fontsize=18,cex = 1.3)),
                                  bottom = textGrob("Time in years",vjust=-2, 
                                                    gp = gpar(fontsize=18,cex = 1.3))) +
  theme(plot.margin=grid::unit(c(0,-10,-10,-10), "mm"))
example.abundance.runawaysim
example.abundance.runawaysim <-  annotate_figure(example.abundance.runawaysim,
                  right = textGrob("Competition", 
                                   check.overlap = T,
                                   vjust = -7.8, hjust = 1.5,
                                   rot=90,
                                   gp = gpar(fontsize=18,cex = 1.3))) +
    theme(plot.margin=grid::unit(c(0,-10,0,0), "mm"))

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



