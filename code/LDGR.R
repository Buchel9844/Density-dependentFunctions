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
##########################################################################################################
# 1. Compute Low- Density Growth Rates for wide range of parameters
##########################################################################################################
source("code/PopProjection_toolbox.R")

# define parameters
nsims <- 1000
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


df.prob.coexist <- NULL
df.prob.long<- NULL
df.list.prob <- NULL
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
    
    df.list <- NULL
    df.list <- grwr(par.dat = p, t.num  = 100,
         function.int)
    df.vert <- df.list[["vert"]]

    df <- df.vert%>% # summary for when time = 1
      dplyr::filter(time==1) %>%
      dplyr::select(invader,grwr,Cgrwc,grwrChesson,grwrStouffer) %>%
      unique() %>%
      melt(id.vars="invader")   %>%
      mutate( vars = paste(invader,variable,sep=".")) %>%
      select( value, vars) %>%
      spread( key=vars, value=value) %>%
      mutate(prob.coex.chess = case_when(Ni.grwrChesson >= 0 &
                                Nj.grwrChesson >= 0  ~ 1, # coexistence
                                Ni.grwrChesson < 0 &
                                Nj.grwrChesson < 0  ~ -1,
                                 TRUE~ 0 ), # competitive exclusion
             prob.coex.id.chess = case_when(Ni.grwrChesson >= 0 & Nj.grwrChesson >= 0  ~ "j_i", # coexistence
                                      Ni.grwrChesson < 0 & Nj.grwrChesson > 0 ~ "j", 
                                      Ni.grwrChesson > 0 & Nj.grwrChesson < 0 ~ "i", # competitive exclusion
                                      Ni.grwrChesson < 0 & Nj.grwrChesson < 0  ~ "-1"),#priority effect
             prob.coex.stouff = case_when(Ni.grwrStouffer >= 0 & Nj.grwrStouffer >= 0  ~ 1, # coexistence
                                   Ni.grwrStouffer < 0 &  Nj.grwrStouffer < 0  ~ -1,
                                   TRUE~ 0 )) %>%
      merge(df.n)
    
    
    df.list.prob[[paste0("int_", i,"_funct_",function.int)]] <-  df.list
    df.prob.long <- bind_rows(df.prob.long,df.vert)
    df.prob.coexist <- bind_rows( df.prob.coexist, bind_cols(df.n, df) )

  }
}

# save .csv

save(df.list.prob,
     file= "results/df.list.prob.Rdata")

write_csv(x = df.prob.coexist, col_names = TRUE, 
          file = paste0("results/df.prob.coexist.csv"))
write_csv(x = df.prob.hory, col_names = TRUE, 
          file = paste0("results/df.prob.hory.csv"))
write_csv(x = df.prob.long , col_names = TRUE, 
          file = paste0("results/df.prob.long.csv.gz"))


#check

if(df$Nj.grwrChesson == log(df.long[which(df.long$invader== "Nj"),]$Nj[2])  -
  log(df.long[which(df.long$invader== "Ni"),]$Nj[2]/df.long[which(df.long$invader== "Ni"),]$Nj[1])){
  
  print("Growth when rare compute correctly as: \n
        (GR of Nj when invading Ni) - (GR of Nj when at eq. and Ni invader)")
}

##########################################################################################################
# 2. Visualise population dynamics
##########################################################################################################

df.prob.long[which(df.prob.long$sim.i == 6 &
                     (df.prob.long$function.int == 1 |
                     df.prob.long$function.int == 4)),] %>%
  gather(Ni, Nj, key="species", value="abundance") %>%
  ggplot(aes(x=time,group=species)) + 
  geom_smooth(aes(y=abundance,color=species,fill=species),alpha=0.2,size=0.5, linetype="dashed") +
  geom_line(aes(y= abundance,color=species),alpha=0.8) +
  theme_bw() + facet_wrap(function.int~invader) + 
  geom_text(aes(label=grwrChesson, y= 7, x=75),color=rep(c("black",rep("NA",times=100)),times=8))+ 
              scale_colour_colorblind() +
  scale_fill_colorblind() +
  labs(title="abundances Ni and Nj over time when invading for function 1 and function 4")

ggsave(last_plot(),
       file = "figures/example.alternative.state.pdf")

# for sim = 2 in above loop
GRWR.list <- list()
sim = 6
t.num= 100
for(function.int in c(1:4)){
  par.dat <- params[[sim]]
N3 <- df.prob.long[which(df.prob.long$sim.i == sim &
                           df.prob.long$function.int == function.int &
                           df.prob.long$invader =="Ni"),c("Ni","Nj")][1,] 
N4 <- df.prob.long[which(df.prob.long$sim.i == sim &
                           df.prob.long$function.int == function.int &
                           df.prob.long$invader =="Nj"),c("Ni","Nj")][1,]

GRWR.list[[function.int]] <- data.frame(time=c(0:t.num),
                                        GRWR.i = #log(Ricker_solution_ODE(state = N4, pars= par.dat, gens=t.num,
                    #function.int)[,"Ni"]) -
  log(Ricker_solution_ODE(state =  N3, pars= par.dat, gens=t.num,function.int)[,"dNi"]),
  GRWR.j = log(Ricker_solution_ODE(state = N4, pars= par.dat, gens=t.num,
                                   function.int)[,"dNj"]) ) %>%
  mutate(coex = case_when(GRWR.i > 0 & GRWR.j> 0 ~ 0,
                          GRWR.i == 0 & GRWR.j == 0 ~ 0),
         col.coex = case_when(GRWR.i >0 & GRWR.j> 0 ~ "coex",
                          GRWR.i == 0 & GRWR.j == 0 ~ "stable")) %>%
  gather(GRWR.i, GRWR.j, key="species", value="GRWR") %>%
  ggplot(aes(x=time)) + 
  geom_line(aes(y=GRWR,color=species),alpha=0.7,size=1)  + 
  geom_point(aes(y=coex,color=col.coex ),shape=1,size= rep(c(1,2,rep(1,each=99)),time=2)) + 
  geom_vline(aes(xintercept=1),linetype="dashed",color="black",alpha=0.7) +
  theme_bw() + scale_colour_colorblind() +
  guides(alpha="none")+
  labs(title=paste0("GRWR at each time step for Ni and Nj \n for sim ",sim," and function ",function.int)) 
}
ggarrange(plotlist = GRWR.list, ncol=2, nrow=2,common.legend = T,
          legend = "right")

ggsave(last_plot(),
       file = "figures/example.GRWR.all.timesteps.pdf")


    
##########################################################################################################
# 3. Compute slopes of growth over time
##########################################################################################################
 
df.glm_all <- NULL

for(i in 1:nsims){
  for( function.int in 1:4){
    print(paste0("int ", i,"for funct ",function.int))
    
df.Ni <- df.prob.long[which(df.prob.long$sim.i == i &
                     df.prob.long$function.int == function.int &
                     df.prob.long$invader =="Ni"),]
if(sum(is.na(df.Ni$Ni))>1) next
if(df.Ni$Nj[1] ==0) next
if(df.Ni$Ni[2] >1000 | df.Ni$Nj[2] >1000) next

df.Ni.glm <- as.data.frame(summary(glm(formula = Ni ~ time, data = df.Ni,
            family = "gaussian"))$coefficients) %>%
  rownames_to_column(var="coeff") %>%
  mutate(coeff = case_when(coeff =="(Intercept)" ~ "intercept",
                           coeff =="time"~"time")) %>%
  select(coeff,Estimate ) %>%
  mutate(coeff = paste0( coeff,".Ni"))%>%
  spread(coeff,Estimate) %>%
  mutate(com.comp.Ni = intercept.Ni + 100*time.Ni,
         com.comp.prob.Ni = case_when(com.comp.Ni>0~T,
                                   com.comp.Ni<0~F),
         sim.i = i,
         function.int = function.int)


df.Nj <- df.prob.long[which(df.prob.long$sim.i == i &
                              df.prob.long$function.int == function.int &
                              df.prob.long$invader =="Nj"),]
if(sum(is.na(df.Nj$Nj))>1) next
if(df.Nj$Ni[1] ==0) next
if(df.Nj$Ni[2] >1000 | df.Nj$Nj[2] >1000) next
df.Nj.glm <- as.data.frame(summary(glm(formula = Nj ~ time, data = df.Nj,
                                       family = "gaussian"))$coefficients) %>%
  rownames_to_column(var="coeff") %>%
  mutate(coeff = case_when(coeff =="(Intercept)" ~ "intercept",
                           coeff =="time"~"time")) %>%
  select(coeff,Estimate ) %>%
  mutate(coeff = paste0( coeff,".Nj"))%>%
  spread(coeff,Estimate) %>%
  mutate(com.comp.Nj = intercept.Nj + 100*time.Nj,
         com.comp.prob.Nj = case_when(com.comp.Nj>0~T,
                                      com.comp.Nj<0~F),
         sim.j = i,
         function.int = function.int)

df.glm <- right_join(df.Ni.glm,df.Nj.glm) %>%
  mutate(com.comp.coex = case_when(com.comp.Ni > 0 &com.comp.Nj > 0 ~ "j_i",
                          com.comp.Ni <= 0 & com.comp.Nj > 0 ~ "j",
                          com.comp.Ni > 0 & com.comp.Nj <= 0 ~ "i",
                          com.comp.Ni <= 0 & com.comp.Nj <= 0 ~ "-1"))

df.glm_all <- bind_rows(df.glm_all,df.glm )
 }
}
head(df.glm_all)
str(df.glm_all)
ggplot(df.glm_all,aes(x=time.Ni, 
                      y=time.Nj,color=as.factor(function.int)))+ geom_point(alpha=0.5) + theme_bw() 




 ##########################################################################################################
# 4. Standardisation
##########################################################################################################


df.prob.coexist <- read.csv("results/df.prob.coexist.csv")
df.prob.coexist.std <- df.prob.coexist %>%
  mutate(prob.coex.chesson = case_when(Ni.grwrChesson > 0 & Nj.grwrChesson > 0  ~ 1, # coexistence
                                 Ni.grwrChesson <= 0 & Nj.grwrChesson > 0 | 
                                 Ni.grwrChesson > 0 & Nj.grwrChesson <= 0 | 
                                 Ni.grwrChesson == 0 & Nj.grwrChesson == 0 ~ 0, # competitive exclusion
                                 Ni.grwrChesson < 0 & Nj.grwrChesson < 0  ~ 0)) %>%
  left_join(df.glm_all, by=c("sim.i","sim.j","function.int")) %>%
  mutate(prob.com.comp.coex =case_when(com.comp.Ni > 0 &com.comp.Nj > 0 ~ 1,
                                       com.comp.Ni <= 0 &com.comp.Nj > 0 ~ 0,
                                       com.comp.Ni > 0 &com.comp.Nj <= 0 ~ 0,
                                       com.comp.Ni <= 0 &com.comp.Nj <= 0 ~ 0))


stand.variable <- c(paste0("Nmax",c(".i.j",".i.i",".j.i",".j.j")),
                    paste0("c",c(".i.j",".i.i",".j.i",".j.j")),
                    paste0("a_initial",c(".i.j",".i.i",".j.i",".j.j")),
                    paste0("a_slope",c(".i.j",".i.i",".j.i",".j.j"))
                    )

df.prob.coexist.std[,stand.variable] <- lapply(df.prob.coexist.std[,stand.variable],
                                scale)
df.prob.coexist.std[,"prob.coex"] <- as.numeric(df.prob.coexist.std[,"prob.coex"])
str(df.prob.coexist.std)

ggplot( df.prob.coexist.std, aes(x=prob.coex.id, 
                                 fill=as.factor(function.int)))+
  stat_count() +theme_bw()
ggplot( df.prob.coexist.std, aes(x=com.comp.coex , 
                             fill=as.factor(function.int)))+
  stat_count() +theme_bw()
       
##########################################################################################################
# 5. Sensitivity analysis
##########################################################################################################

# Incorporate interaction terms

model.0 <- glm(prob.com.comp.coex ~ as.factor(function.int), 
               data = df.prob.coexist.std[which(df.prob.coexist.std$prob.coex >= 0),],
               family = "binomial")

model.1 <- glm(formula(paste0("prob.com.comp.coex  ~ ",
                             paste(paste0("a_initial",c(".i.j",".i.i",".j.i",".j.j")),collapse = "+"))), 
   data = subset(df.prob.coexist.std,function.int==1 & prob.coex >= 0),
   family = "binomial")

model.2 <- glm(formula(paste0("prob.com.comp.coex  ~ ",
                             paste(c(paste0("Nmax",c(".i.j",".i.i",".j.i",".j.j")),
                                     paste0("a_initial",c(".i.j",".i.i",".j.i",".j.j")),
                                     paste0("a_slope",c(".i.j",".i.i",".j.i",".j.j"))),collapse = "+"))), 
              data = subset(df.prob.coexist.std,function.int==2 & prob.coex >= 0),
              family = "binomial")

model.3 <- glm(formula(paste0("prob.com.comp.coex  ~ ",
                             paste(stand.variable,collapse = "+"))), 
              data = subset(df.prob.coexist.std,function.int==3 & prob.coex >= 0),
              family = "binomial")

model.4 <- glm(formula(paste0("prob.com.comp.coex  ~ ",
                             paste(stand.variable,collapse = "+"))), 
              data = subset(df.prob.coexist.std,function.int==4 & prob.coex >= 0),
              family = "binomial")

# build dataframe to plot
ggplot( df.prob.coexist.std, aes(x=prob.coex.id, 
                             fill=as.factor(function.int)))+
  stat_count() +theme_bw()


sens_out <- tidy(model.1) %>% 
  mutate(function.int = 1) %>% 
  bind_rows(
    (tidy(model.0) %>%
       mutate(term = "function",
              function.int = c(1,2,3,4))),
    (tidy(model.2) %>%
      mutate(function.int = 2)),
    (tidy(model.3) %>%
  mutate(function.int = 3)),
    (tidy(model.4) %>%
      mutate(function.int = 4))
    ) %>%
  filter(term != "(Intercept)")

sens_out$term <- factor(sens_out$term, 
                        levels = 
                          c("function",
                            paste0("a_initial",c(".i.i",".j.j",".i.j",".j.i")),
                            paste0("a_slope",c(".i.i",".j.j",".i.j",".j.i")),
                            paste0("Nmax",c(".i.i",".j.j",".i.j",".j.i")),
                            paste0("c",c(".i.i",".j.j",".i.j",".j.i"))
                          )
)
# change to have the focal species first
x_axis_labs <- c("function",
                 paste0("a_initial",c(".i.i",".j.j",".j.i",".i.j")),
                 paste0("a_slope",c(".i.i",".j.j",".j.i",".i.j")),
                 paste0("Nmax",c(".i.i",".j.j",".j.i",".i.j")),
                 paste0("c",c(".i.i",".j.j",".j.i",".i.j"))
                 )

my_cols <- qualitative_hcl(n =4)

int_sensitivity_plot <- sens_out %>% 
  
  ggplot(aes(y = estimate, ymin = (estimate - std.error), ymax = (estimate + std.error),
             color = as.factor(function.int), 
             fill = as.factor(function.int), x = term)) +
  #geom_vline(xintercept = c(1:9+0.5), alpha = 0.5, color = "gray50", size = 0.2) +
  geom_bar(stat = "identity",
           position = position_dodge(), alpha = 0.7) +
  geom_errorbar(position = position_dodge(.9), width = 0.5, show.legend = FALSE) +
  scale_color_manual(values = darken(my_cols, amount = .2)) + 
  scale_fill_manual(values = my_cols) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = c(.07,.92),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.key.size = unit(0.1, units = "in"),
        axis.text.x = element_text(size = 12, angle = 45, 
                                   hjust = 1, margin = margin(t =10))) +
  #facet_wrap(~comp, nrow = 2) +
  labs(x = "", y = "Effect Size on \n probability of coexistence", fill = "function", color = "function") +
  scale_x_discrete(labels = x_axis_labs)

int_sensitivity_plot
ggsave("figures/sensitivity_analysis.pdf", width = 10, height = 8)



# plot probability of each function to have coexistence

df <- df.prob.coexist.std

fits <- predict(model.0, newdata = df, se.fit = TRUE)

# Get odds from modrl
df$prediction <- exp(fits$fit) 
df$upper <- exp(fits$fit + 1.96 * fits$se.fit)
df$lower <- exp(fits$fit - 1.96 * fits$se.fit)

# Convert odds to probabilities
df$prediction <- df$prediction / (1 + df$prediction)
df$upper <- df$upper / (1 + df$upper)
df$lower <- df$lower / (1 + df$lower)

# Plot probabilities

ggplot(df, aes(function.int, prediction)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.25, size = 1, position = position_dodge(width = 0.4)) +
  geom_point(shape = 21, size = 3, 
             position = position_dodge(width = 0.4)) +
  theme_light(base_size = 16) +
  scale_y_continuous(name = "Probability of coexistence", limits = c(0, 1),
                     labels = scales::percent)

