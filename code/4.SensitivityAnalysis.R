# Script for sensitivity analysis
library(forecast) # for ARIMA function
library(lmtest) # for ARIMA eigen value
library(ppcor) # for variable correlation
library(tsvr) # for synchrony 
library(tidyverse) # for synchrony 
library(ggplot2) 
library(ggthemes)
library(ggpattern)
library(wesanderson)
library(colorspace)
library(broom)
library(truncnorm)
###########################################################################################################
# 1. Compute Low- Density Growth Rates for wide range of parameters without coexistence restriction
##########################################################################################################
source("code/PopProjection_toolbox.R")

# define parameters
nsims <- 2000
species <- 2
function.int <- c(1:4)
t.num = 100
gen = seq(from=0, to=t.num , by=1)
# run sim
set.seed(1608)

# Set parameters lists

white.noise <- function(t.num){
  e_noise = rnorm(t.num+1,mean=0, sd=0.1) 
  for( n in 1:sample(2:10)[1]){
    e_noise[n] <- rnorm(1,mean=0,
                        sd=abs(rnorm(1, mean=0, sd=1)))
  }
  return( e_noise)
}

params.sens <-list()
for( i in 1:nsims){
  params.sens[[i]] <- list(
    # stable parameters
    sim= i, 
    g = c(gi=0.9,gj=0.9), # germination rate of seeds
    s = c(si= 0.9,sj= 0.9), # survival rate of seeds
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
                       data = runif(n=species*species, min=-1, max=0),
                       dimnames = list(c("i", "j"),
                                       c("i", "j"))),
    
    c =matrix(ncol=species, nrow = species, # a_slope, only positive values
              data =runif(n=species*species, min=-1, max=0),
              dimnames = list(c("i", "j"),
                              c("i", "j"))),
    
    e_seasonal = sin((2*pi/20)*gen)*abs(rnorm(1,mean=0, sd=0.1)),
    e_noise = white.noise(t.num)
    
  )
}


df.sim.sens  <- NULL

t.num = 100 # number of generation
set.seed(1608)
for(i in 1:nsims){
  for( function.int in 1:4){
    for(add_external_factor in c("none")){
      print(paste0("int ", i,"for funct ",function.int, add_external_factor))
      df.n <- NULL
      p <- params.sens[[i]]
      df.n <- as.data.frame(params.sens[[i]][c("sim","Nmax","a_slope","a_initial","c")]) %>%
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
      
      df.initc <- Ricker_solution(state = c(1,1), pars= p, gens=t.num,
                                  function.int,add_external_factor) %>% # simulation of initial condition = 1,1
        mutate(invader = "both") %>%
        bind_cols(df.n) %>% 
        bind_rows(df.inv)%>% 
        mutate(external_factor= add_external_factor)
      
      for(inv in c("both","Ni","Nj")){ # compute stability metric
        df.sens.n <- df.initc %>%
          filter(invader == inv) %>%
          mutate( stability.i = mean(Ni)/var(Ni),
                  stability.j = mean(Nj)/var(Nj)) %>%
          mutate( stability.i = case_when(mean(Ni) < 0.10~0,
                                          T~stability.i),
                  stability.j = case_when(mean(Nj) < 0.10~0,
                                          T~stability.j))
        
        df.sim.sens <- bind_rows(df.sim.sens, df.sens.n)
      }
    }
  }
}

# save .csv
df.sim.sens <- df.sim.sens %>%
  mutate(function.name = case_when(function.int==1 ~"1.Constant",
                                   function.int==2 ~"2.Linear",
                                   function.int==3 ~"3.Exp",
                                   function.int==4 ~"4.Sigmoid"
                                   ))

write.csv(df.sim.sens , 
          file = paste0("results/df.sim_sens.csv.gz"))
load("results/df.sim.csv.gz")
df.sim.sens = fread("results/df.sim_sen.csv.gz")
head(df.sim)

###########################################################################################################
# 3. Visualise of abundance over time with stability
##########################################################################################################
sim = 100
t.num= 100
#---- Abundance through time with external factors----

example.abundance.sens <- df.sim.sens[which(df.sim.sens$sim.i == sim &
                                                       df.sim.sens$invader == "both"),] %>%
  gather(Ni, Nj, key=species, value=abundance) %>%
  ggplot(aes(x=time)) +
  geom_line(aes(y= abundance,color=species)) +
  geom_line(aes(y= external.fact),color="black") +
  theme_bw() + 
  scale_color_manual(values=c("#332288", "#999933")) +
  scale_fill_manual(values=c("#332288", "#999933")) +
  facet_wrap(. ~ function.name,nrow=2,ncol=2,
             scales="free") #+
  #annotate(geom="text",aes(label=stability.i),size=1,x=1,y=1) +
  
example.abundance.sens
ggsave(example.abundance.sens,
       file = "figures/example.abundance.sens.pdf")
###########################################################################################################
# 3. Scale parameters for sensitivity analysis
###########################################################################################################

df.sim.sens.std.unscale <- df.sim.sens %>%
  filter(time==100)
df.sim.sens.std.unscale <- data.frame(function.int = c(df.sim.sens.std.unscale$function.int,df.sim.sens.std.unscale$function.int),
           function.name =c(df.sim.sens.std.unscale$function.name,df.sim.sens.std.unscale$function.name),
           stability =c(df.sim.sens.std.unscale$stability.i,df.sim.sens.std.unscale$stability.j),
           alpha.decay.intra = c(df.sim.sens.std.unscale$a_slope.i.i,df.sim.sens.std.unscale$a_slope.j.j),
           alpha.decay.inter = c(df.sim.sens.std.unscale$a_slope.i.j,df.sim.sens.std.unscale$a_slope.j.i),
           alpha.0.intra = c(df.sim.sens.std.unscale$a_initial.i.i,df.sim.sens.std.unscale$a_initial.j.j),
           alpha.0.inter = c(df.sim.sens.std.unscale$a_initial.i.j,df.sim.sens.std.unscale$a_initial.j.i),
           C.intra = c(df.sim.sens.std.unscale$c.i.i,df.sim.sens.std.unscale$c.j.j), 
           C.inter = c(df.sim.sens.std.unscale$c.i.j,df.sim.sens.std.unscale$c.j.i),
           N.opt.intra = c(df.sim.sens.std.unscale$Nmax.i.i,df.sim.sens.std.unscale$Nmax.j.j),
           N.opt.inter = c(df.sim.sens.std.unscale$Nmax.i.j,df.sim.sens.std.unscale$Nmax.j.i))

ggplot(df.sim.sens.std.unscale, aes(x=stability, group=function.int, color=function.int ))+ geom_density()
stand.variable <- c("alpha.decay.intra","alpha.decay.inter",
                    "alpha.0.intra","alpha.0.inter",
                    "C.intra","C.inter",
                    "N.opt.intra","N.opt.inter")
names(df.sim.sens.std.unscale)
df.sim.sens.std <- df.sim.sens.std.unscale
df.sim.sens.std[,stand.variable] <- lapply(df.sim.sens.std.unscale[,stand.variable],
                                      scale)

###########################################################################################################
# 4. Sensitivity analysis
##########################################################################################################
#---- Sensitivity analysis -----
df.sim.sens.std <- df.sim.sens.std %>%
  mutate(stability = case_when(stability > 100 ~ 100,
                               T~stability),
         function.int = as.factor(function.int))

model.0 <- lm(stability ~ as.factor(function.int),
              df.sim.sens.std)

model.1 <- glm(formula(paste0("stability ~ ",
                             paste(c("alpha.0.intra","alpha.0.inter"),collapse = "+"))), 
              data = df.sim.sens.std[which(df.sim.sens.std$function.int==1),])


model.2 <- glm(formula(paste0("stability~ ",
                             paste(c("alpha.0.intra","alpha.0.inter",
                                     "alpha.decay.intra","alpha.decay.inter"),collapse = "+"))), 
              data = df.sim.sens.std[which(df.sim.sens.std$function.int==2),])

model.3 <- glm(formula(paste0("stability ~ ",
                             paste(c("alpha.0.intra","alpha.0.inter",
                                     "alpha.decay.intra","alpha.decay.inter",
                                     "C.intra","C.inter",
                                     "N.opt.intra","N.opt.inter"),collapse = "+"))), 
              data = df.sim.sens.std[which(df.sim.sens.std$function.int==3),])

model.4 <- glm(formula(paste0("stability ~ ",
                             paste(c("alpha.0.intra","alpha.0.inter",
                                     "alpha.decay.intra","alpha.decay.inter",
                                     "C.intra","C.inter",
                                     "N.opt.intra","N.opt.inter"),collapse = "+"))), 
              data = df.sim.sens.std[which(df.sim.sens.std$function.int==4),])

predict.1 <- bind_cols(df.sim.sens.std.unscale[which(df.sim.sens.std.unscale$function.int==1),],
            data.frame(predict.stability=predict(model.1,df.sim.sens.std[which(df.sim.sens.std$function.int==1),]))) %>%
  gather(c("alpha.0.intra","alpha.0.inter"),key="term",value="value")



predict.2 <- bind_cols(df.sim.sens.std.unscale[which(df.sim.sens.std.unscale$function.int==2),],
                       data.frame(predict.stability=predict(model.2,df.sim.sens.std[which(df.sim.sens.std$function.int==2),]))) %>%
  gather(c("alpha.0.intra","alpha.0.inter",
           "alpha.decay.intra","alpha.decay.inter"),key="term",value="value")

predict.3 <- bind_cols(df.sim.sens.std.unscale[which(df.sim.sens.std.unscale$function.int==3),],
                       data.frame(predict.stability=predict(model.3,df.sim.sens.std[which(df.sim.sens.std$function.int==3),]))) %>%
  gather(c("alpha.0.intra","alpha.0.inter",
           "alpha.decay.intra","alpha.decay.inter",
           "C.intra","C.inter",
           "N.opt.intra","N.opt.inter"),key="term",value="value")

predict.4 <- bind_cols(df.sim.sens.std.unscale[which(df.sim.sens.std.unscale$function.int==4),],
                       data.frame(predict.stability=predict(model.4,df.sim.sens.std[which(df.sim.sens.std$function.int==4),]))) %>%
  gather(c("alpha.0.intra","alpha.0.inter",
           "alpha.decay.intra","alpha.decay.inter",
           "C.intra","C.inter",
           "N.opt.intra","N.opt.inter"),key="term",value="value")

predict <- bind_rows(predict.1,predict.2,predict.3,predict.4)

my_cols <- c("#AA4499", "#DDCC77","#88CCEE", "#44AA99")


predict_sensitivity_plot <- predict %>%
  ggplot(aes(x=value,color=function.name,fill=function.name)) +
  #geom_smooth(aes(y=predict.stability), se=FALSE)+
  geom_point(aes(y=stability),alpha=0.1,size=0.5) +
  geom_smooth(aes(y=predict.stability))+
  theme_bw()+
  scale_color_manual(values = darken(my_cols, amount = .1)) + 
  scale_fill_manual(values = my_cols) + 
  facet_wrap(~term,nrow = 2,dir="v",scales="free") +
  scale_y_continuous(limits=c(0,100))+
  theme(legend.background = element_rect(color = "white"),
        legend.position ="right",
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, units = "in"))+
  #facet_wrap(~comp, nrow = 2) +
  labs(title="Prediction of parameters' effect for each function on the stability of one species in 2-species community",
       y = "Stability ratio \nlow values equal lower stability",
       x = "Range of value specific to the parameters, from low to high",
       fill = "function", color = "function") 

predict_sensitivity_plot
ggsave("figures/predict_sensitivity_plot.pdf", 
       width = 10, height = 8)


sens_out <- tidy(model.4) %>% 
  mutate(function.name = "4.Sigmoid") %>% 
  bind_rows(
    #(tidy(model.0) %>%
     #  mutate(term = "function",
     #         function.int = c(1,2,3,4))),
    (tidy(model.1) %>%
       mutate(function.name = "1.Constant")),
    (tidy(model.3) %>%
       mutate(function.name = "3.Exp")),
    (tidy(model.2) %>%
       mutate(function.name = "2.Linear"))
  ) %>%
  filter(term != "(Intercept)")

sens_out$term <- factor(sens_out$term, 
                        levels = 
                          c("alpha.0.intra","alpha.0.inter",
                            "alpha.decay.intra","alpha.decay.inter",
                            "C.intra","C.inter",
                            "N.opt.intra","N.opt.inter")
)
        


int_sensitivity_plot <- sens_out %>% 
  filter(term!="function") %>%
  ggplot(aes(x = estimate, xmin = (estimate - std.error), xmax = (estimate + std.error),
             color = as.factor(function.name), 
             fill = as.factor(function.name), y = as.factor(term))) +
  #geom_vline(xintercept = c(1:9+0.5), alpha = 0.5, color = "gray50", size = 0.2) +
  geom_bar(stat = "identity",width =0.6,
           position = position_dodge(), alpha = 0.7) +
  geom_errorbar(position = position_dodge(), width = 0.6, show.legend = FALSE) +
  scale_color_manual(values = darken(my_cols, amount = .2)) + 
  scale_fill_manual(values = my_cols) +
  #scale_x_continuous(labes=("lowe"))
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = c(.87,.90),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(0.1, units = "in"),
        axis.text.y = element_text(size = 12, angle = 0, 
                                   hjust = 1, margin = margin(t =10))) +
  #facet_wrap(~comp, nrow = 2) +
  labs(title="Effect size of parameters for each function on the stability of one species in 2-species community",
       y = "", x = "Effect Size on Stability of i \n ratio mean/var", fill = "", color = "") 

int_sensitivity_plot
ggsave("figures/sensitivity_analysis_stability.pdf", 
       width = 10, height = 8)

# plot probability of each function to have coexistence

df <- df.sim.sens.std

fits <- predict(model.0, newdata = df, se.fit = TRUE)

# Get odds from modrl
df$prediction <- fits$fit
df$upper <- fits$fit + 1.96 * fits$se.fit
df$lower <- fits$fit - 1.96 * fits$se.fit

# Plot probabilities

prediction.coexistence <- ggplot(df, aes(x=as.factor(function.name), 
                                         y=prediction,
                                         color=as.factor(function.name))) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.25, size = 1, position = position_dodge(width = 0.4)) +
  geom_point(shape = 21, size = 3, 
             position = position_dodge(width = 0.4)) +
  geom_hline(yintercept=0,color="black",linetype="dashed")+
  scale_color_manual(values=my_cols) +
  guides(color="none")+
  theme_light(base_size = 16) +
  labs(color="",x="",title="Stability distribution for one species in a 2-species community,\nfrom unstable (neg) to stable (positive)",
       y="Stability prediction")
prediction.coexistence
ggsave(prediction.coexistence,
       "figures/prediction.stability.pdf")
