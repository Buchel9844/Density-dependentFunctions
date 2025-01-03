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
##########################################################################################################
# Stability metric
##########################################################################################################
 source("code/3.1.TimeSerie_toolbox.R")
df.stability <- as.data.frame(df.sim.comcomp)
 summary.df.stability <- NULL
 nsims <- 1500
 t.num = 200 # number of generation
 stand.variable <- c(paste0("a_initial",c(".i.i",".j.j",".i.j",".j.i")),
                     paste0("a_slope",c(".i.i",".j.j",".i.j",".j.i")),
                     paste0("c",c(".i.i",".j.j",".i.j",".j.i")),
                     paste0("Nmax",c(".i.i",".j.j",".i.j",".j.i")))
 
 
 df.stability.summary <- NULL
 for(i in 1:nsims){
   for( function.int in 1:4){
     for(add_external_factor in c("No external factor","Noisy change","Periodic change")){
       print(paste0("int ", i,"for funct ",function.int, add_external_factor))
       
   df.stability.n <-  as.data.frame(df.stability[which(df.stability$sim.i == i & 
                                           df.stability$time > 99 & 
                                           df.stability$function.int == function.int &
                                           df.stability$external_fact == add_external_factor),] )
   comp.com.n <- df.stability.n$comp.com[1]
   if(!comp.com.n  %in% c("one-species community","two-species community")) next
#Mean  over time   
   mean.i=c(mean(df.stability.n$Ni))
   mean.j=c(mean(df.stability.n$Nj))
   
 #Mean Variance over time
 # Plaza et al,2012 https://doi.org/10.1111/j.1654-1103.2011.01381.x
 #  amplitude of population fluctuations by means of the standard deviation
   # ratio between median and interquantile
   stability.int <- c()
   if(exists("stability.significance")){rm("stability.significance")}
   
   ratio_i <- mean(df.stability.n$Ni,na.rm=T)/var(df.stability.n$Ni,na.rm=T) 
   ratio_j <- mean(df.stability.n$Nj,na.rm=T)/var(df.stability.n$Nj,na.rm=T)
   
      msdi <- sd(log10(df.stability.n$Ni))
      msdj <- sd(log10(df.stability.n$Nj))
      # for j
    if(is.na(ratio_j)){
      stability.int <- 0
    }
    if(!is.na(ratio_j)){
    if(is.infinite(ratio_j)|ratio_j > 1){
      stability.int <- 1
    }else{
      stability.int <- 0
    }
    }
      #for i
      if(is.na(ratio_i)){
        stability.int <- stability.int + 0
      }
      if(!is.na(ratio_i)){
        if(is.infinite(ratio_i)|ratio_i > 1){
        stability.int <- stability.int + 1
        }else{
          stability.int <- 0
        }
      }
      
      if(stability.int ==0){
        stability.significance = "non-equilibrium"
      }
      if(stability.int ==1){
        stability.significance = "one-species equilibrium"
      }
      if(stability.int ==2|(stability.int ==1 & comp.com.n=="one-species community")){
        stability.significance = "stable"
      }
      
      if(!exists("stability.significance")){
        stability.significance = "no"
        ratio_i <- NA
        ratio_j <- NA
        msdi <- NA
        msdj <- NA
       }
#ARIMA is the abbreviation for AutoRegressive Integrated Moving Average.      
      #for j
      if(exists(c("p.j"))|exists(c("p.i"))){rm(p.j,q.j,p.i,q.i)}
   
      if(median(df.stability.n$Nj) > 0 & var(df.stability.n$Nj) > 0  & !is.na(mean(log(df.stability.n$Nj)))){
      fitARIMA.j <- auto.arima(df.stability.n$Nj) # function in R uses a combination of unit root tests, minimization of the AIC and MLE to obtain an ARIMA model
      p.j = summary(fitARIMA.j)$arma[1]#ARIMA(p,d,q)
      q.j = summary(fitARIMA.j)$arma[2]
      # summary(fitARIMA.j)
      if(p.j > 0 & q.j >= 0){
       coeff.j = summary(fitARIMA.j)$coef
       arima.extract.j  <- ARMApqREMLfunct(coeff.j ,log(df.stability.n$Nj),p.j,q.j)
       names(coeff.j) <- paste0(names(coeff.j), ".j")
       coeff.j <- as.data.frame(t(as.matrix(coeff.j)))
       coeff.j$p.j <- as.integer(p.j) #ARIMA(p,d,q)
       coeff.j$q.j = q.j
       coeff.j$eigB.j <- arima.extract.j$eigB
       coeff.j$LL.j <- arima.extract.j$LL
     }else{
       coeff.j <- data.frame(p.j = p.j, q.j = q.j,
                           eigB.j = 0)
     }
      }else{
        coeff.j <- data.frame(p.j = NA, q.j = NA,
                              eigB.j = NA)
      }
      # for i 
        if(median(df.stability.n$Ni) > 0.001 & var(df.stability.n$Ni) > 0.001 & !is.na(mean(log(df.stability.n$Ni)))){
      fitARIMA.i <- auto.arima(df.stability.n$Ni)
      p.i = summary(fitARIMA.i)$arma[1]#ARIMA(p,d,q)
      q.i = summary(fitARIMA.i)$arma[2]
      
      if(p.i > 0 & q.i >= 0){
        coeff.i = summary(fitARIMA.i)$coef
        arima.extract.i  <- ARMApqREMLfunct(coeff.i, log(df.stability.n$Ni),p.i,q.i)
        names(coeff.i) <- paste0(names(coeff.i), ".i")
        coeff.i <- as.data.frame(t(as.matrix(coeff.i)))
        coeff.i$p.i = p.i#ARIMA(p,d,q)
        coeff.i$q.i = q.i
        coeff.i$eigB.i <- arima.extract.i$eigB
        coeff.i$LL.i <- arima.extract.i$LL
      }else{
        coeff.i <- data.frame(p.i = p.i, q.i = q.i,
                              eigB.i = 0)
        }
        }else{
          coeff.i <- data.frame(p.i = NA, q.i = NA,
                                eigB.i = NA)
        }
      
      if(exists("oscillation.significance")){rm("oscillation.significance")}
      if( (!exists(c("p.j")) | !exists(c("p.i")))){
        oscillation.significance = "no"
      }else{ 
        if( ((coeff.j$p.j==0|coeff.j$q.j==0|
              is.na(coeff.j$p.j)|is.na(coeff.j$q.j)) & 
              (coeff.i$p.i > 0|coeff.i$q.i>0) & 
              coeff.i$eigB.i > 0 
              ) |
            ( (coeff.i$p.i==0|coeff.i$q.i==0 |
              is.na(coeff.i$p.i)|is.na(coeff.i$q.i)) &
              (coeff.j$p.j >0|coeff.j$q.j>0 ) &
              coeff.j$eigB.j > 0 )
            ){
          oscillation.significance = "half_oscillation"
          if(comp.com.n=="one-species community"){
            oscillation.significance = "oscillation" 
          }
          
        }else{
                if((coeff.i$eigB.i== 0 | coeff.j$eigB.j == 0 |is.na(coeff.i$eigB.i)| is.na(coeff.j$eigB.j)) &
                 (p.j ==0|q.j==0 | is.na(p.j)|is.na(q.j) ) &
                 (p.i ==0|q.i==0 | is.na(p.i)|is.na(q.i))){
                oscillation.significance = "no"
              }else{
                oscillation.significance = "oscillation"
              }
          }
        }
         
        # Multiple conditions with and
 
        # Output
        #[1] "a value is in between 40 and 60"
     
#SYNCHRONIE OF NI AND NJ - > 1 synchrony | < 1 asynchrony
   if(exists("vr.trial")){rm("vr.trial")}
   
     if(   msdi > 0.05  & msdj > 0.05 &!is.na(msdi) &!is.na(msdj) &
           mean(df.stability.n$Ni[c(10:nrow( df.stability.n))]) > 0.05  & mean(df.stability.n$Nj) > 0.05 ){ # otherwise thre are no variation in the data
     df.stability.n.small <- t(as.matrix(data.frame(Ni= df.stability.n$Ni,
                                                    Nj =df.stability.n$Nj)[c(10:nrow( df.stability.n)),]))
     an.error.occured <- FALSE
     tryCatch( { vr.trial <- tsvreq_classic(df.stability.n.small); print(res) }
               , error = function(e) {an.error.occured <<- TRUE})
     if(exists("vr.trial")){
     aggresShort <- aggts(vr.trial, vr.trial$ts[vr.trial$ts<4])[[3]]
     aggresLong <- aggts(vr.trial, vr.trial$ts[vr.trial$ts>=4])[[3]]
     if(aggresShort>1|aggresLong> 1){
       synchrony.significance <- "Long and short synchrony"
       if(aggresShort>1 & aggresLong < 1){
         synchrony.significance <-"Short synchrony"
       }
       if(aggresShort<1 & aggresLong > 1){
         synchrony.significance <-"Long synchrony"
       }
     }else{synchrony.significance <- "No"}
     }else{
       aggresLong <- NA
       aggresShort <- NA
       synchrony.significance <- NA
     }
     }else{
       aggresLong <- NA
       aggresShort <- NA
       synchrony.significance <- NA
     }
   
     df.stability.summary.n <- data.frame(sim = as.integer(i),
                                     function.int = as.integer(function.int),
                                     external_factor = add_external_factor,
                                     mean.i=c(mean(df.stability.n$Ni)),
                                     mean.j=c(mean(df.stability.n$Nj)),
                                     median.i=c(median(df.stability.n$Ni)),
                                     median.j=c(median(df.stability.n$Nj)),
                                     ratio.i = c(ratio_i),
                                     ratio.j = c(ratio_j),
                                     oscillation.significance =   oscillation.significance, 
                                     stability.significance = stability.significance,
                                     synchrony.significance = synchrony.significance,
                                     synchrony.long = aggresLong,
                                     synchrony.short = aggresShort) 

     df.stability.summary.n <- bind_cols(df.stability.summary.n,bind_cols(coeff.j,coeff.i))
     df.stability.summary.n <- full_join(df.stability.summary.n , 
                                         as.data.frame(df.stability.n[1,]), 
                 by = c("sim", "function.int", "external_factor"))
    
    df.stability.summary <- bind_rows(df.stability.summary,df.stability.summary.n)
                                     
 
     }
   }
}

save(file="results/df.stability.summary.csv.gz",
     df.stability.summary) 
 
load("results/df.stability.summary.csv.gz") 
#---- global visualisation -----

df.stability.summary.small <- df.stability.summary

df.stability.summary.small$significance <- apply(df.stability.summary.small[, c("synchrony.significance",
                                              "oscillation.significance")],
                         1, function(x) toString(na.omit(x)))

df.stability.summary.small <- df.stability.summary.small %>%
  mutate(significance = case_when(   oscillation.significance == "half_oscillation" &
                                       synchrony.significance == "no" ~ "one species oscillation",
                                     oscillation.significance == "oscillation" &
                                       synchrony.significance == "no" ~ "oscillation",
                                     synchrony.significance == "Long and short synchrony"&
                                       oscillation.significance == "no"  ~ "Long and short synchrony",
                                     synchrony.significance == "Short synchrony"&
                                       oscillation.significance == "no"  ~ "Short synchrony",
                                     synchrony.significance == "Long synchrony" &
                                       oscillation.significance == "no" ~ "Long synchrony",
                                     ((synchrony.significance =="Long synchrony"|
                                        synchrony.significance =="Short synchrony"|
                                        synchrony.significance =="Long and short synchrony") & 
                                       (oscillation.significance == "oscillation" |
                                          oscillation.significance == "half_oscillation"))  ~ "Synchronous Oscillation",
                                     TRUE ~ "no detected dynamics"
                                  ),
         stability.significance = case_when(stability.significance == "stable" ~"All species",
                                            stability.significance == "one-species equilibrium" ~"1 species out of 2",
                                            stability.significance == "non-equilibrium" ~"0 species")) %>%
  mutate(significance = factor(significance, 
                               levels=c("no detected dynamics","Short synchrony","Long synchrony","Long and short synchrony","Synchronous Oscillation")))

df.stability.summary.small  %>% count(function.int, significance)

levels(as.factor(df.stability.summary.small$significance))
levels(as.factor(df.stability.summary.small$stability.significance))
levels(as.factor(df.stability$external_factor))
levels(as.factor(df.stability$function.name))
df.stability.summary.small[which(is.na(df.stability.summary.small$significance)),]

my_cols <- c("#AA4499", "#DDCC77","#88CCEE", "#44AA99")
cols_interaction <- c("#661100", "#888888", "#6699CC", "#332288") # "#DDCC77","#661100","#117733")
  
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

safe_colorblind_palette_4 <-c("lightgrey","#6699CC","#332288","#117733","#DDCC77")

#colorblind_palette_8 <- safe_colorblind_palette[!safe_colorblind_palette %in% my_cols]
#scales::show_col(safe_colorblind_palette)

df.stability.summary.small <- df.stability.summary.small  %>%
  mutate(function.name = case_when(function.int==1 ~"1.Traditional constant",
                                   function.int==2 ~"2.Linear",
                                   function.int==3 ~"3.Exp",
                                   function.int==4 ~"4.Sigmoid"))

levels(as.factor(df.stability$function.int))
levels(as.factor(df.stability.summary$comp.com))

#---- facet grid ----
df.coexistence.prob <-  df.stability %>%
  filter(external_factor =="No external factor") %>%
  mutate(function.name = case_when(function.int==1 ~"Traditional constant",
                                   function.int==2 ~"Linear",
                                   function.int==3 ~"Exp",
                                   function.int==4 ~"Sigmoid")) %>%
  mutate(function.name  = factor(function.name,
                                 levels=c("Traditional constant","Linear","Exp","Sigmoid"))) %>%
  aggregate(sim ~ function.name + comp.com + Interspecific.interaction, length) %>%
  mutate(sim = round((sim/(1500*100*2))*100,digits=2)) %>%
  spread(Interspecific.interaction,sim) %>%
  replace(is.na(.),0) %>%
  mutate(Total = rowSums(.[,levels(as.factor(df.stability.summary.small$Interspecific.interaction))],na.rm=T))
df.coexistence.prob
names(df.coexistence.prob) <- c("functional form","Community dynamics",
                                   "Competition","One competing,one facilitating",
                                   "Facilitation","Total")
write.csv(df.coexistence.prob ,
          "results/df.coexistence.prob.csv")

df.stability.prob <-  df.stability.summary.small %>%
  filter(external_factor =="No external factor") %>%
  filter(comp.com == "two-species community") %>%
  mutate(function.name = case_when(function.int==1 ~"Traditional constant",
                                   function.int==2 ~"Linear",
                                   function.int==3 ~"Exp",
                                   function.int==4 ~"Sigmoid")) %>%
  mutate(function.name  = factor(function.name,
                                 levels=c("Traditional constant","Linear","Exp","Sigmoid"))) %>%
  mutate(stab.prob = case_when(stability.significance == "All species" & comp.com =="two-species community"~"2 stable species",
                               stability.significance == "All species" & comp.com =="one-species community"~"1 stable species",
                               stability.significance == "1 species out of 2" ~"1 stable species",
                               stability.significance == "0 species" ~"0 species")) %>%
  aggregate(sim ~ function.name + stab.prob + significance,length) %>%
  mutate(sim = round((sim/1500)*100,digits=2))%>%
  spread(significance,sim) %>%
  replace(is.na(.),0) %>%
  mutate(Total = rowSums(.[,levels(as.factor(df.stability.summary.small$significance))],na.rm=T))
df.stability.prob
names(df.stability.prob)[1:2] <- c("functional form","Stability")
write.csv(df.stability.prob,
          "results/df.stability.prob.csv")
plot.stability.prob <- list()
for(fi in 1:4){
  title.vec <- c("Traditional constant",
                 "Linear","Exponential",
                 "Sigmoid")
  y.vec <- c("probability of detection in 2 species-communities",
             "","","")
  y.position <- c(2,0,0,0)
plot.stability.prob[[fi]] <-  df.stability.summary.small %>%
  mutate(function.name = case_when(function.int==1 ~"Traditional constant",
                                   function.int==2 ~"Linear",
                                   function.int==3 ~"Exponential",
                                   function.int==4 ~"Sigmoid")) %>%
  mutate(function.name  = factor(function.name,
                                 levels=c("Traditional constant","Linear",
                                          "Exponential","Sigmoid"))) %>%
  filter(external_factor =="No external factor") %>%
  filter(comp.com == "two-species community") %>%
  mutate(stab.prob = case_when(stability.significance == "All species" & comp.com =="two-species community"~"2 stable species",
                               stability.significance == "All species" & comp.com =="one-species community"~"1 stable species",
                               stability.significance == "1 species out of 2" ~"1 stable species",
                               stability.significance == "0 species" ~"0 species")) %>%
  aggregate(sim ~ stab.prob + significance + function.int+function.name, length) %>%
  mutate(sim = (sim/1500)) %>%
  filter(function.int==fi) %>%
  ggplot() +
  geom_bar(aes(x=as.factor(stab.prob), y=sim, fill=significance ), color="black",
           stat="identity") + 
  #facet_wrap(.~function.name, nrow=1, strip.position= "top") +
  scale_y_continuous(y.vec[fi],
                     labels = scales::percent_format(accuracy = 1),
                     limits=c(0,0.25)) +
  scale_x_discrete("",
                   labels=c("no species\nstable","one species\nstable\n","two species\nstable\n"))+
  scale_fill_manual("Community\ndynamics",
                    labels=c("no detected\ndynamics","short\nsynchrony","long\nsynchrony",
                            "long & short\nsynchrony","sychrony &\noscilattion"),
                    values=  safe_colorblind_palette_4 ) +
  guides(fill= guide_legend(override.aes = list(pattern = "none"),
                           # position =c(.95, .95),
                            direction="horizontal",
                            byrow = TRUE,
                            nrow = 1,
                            title.hjust = 0.1)) +
  labs( title=title.vec[fi]) +
  theme_minimal() +
  theme(panel.spacing.x = unit(10, "mm"),
        legend.position = c(0.5, -0.25),
        panel.grid.major.x  = element_blank(),
        legend.text = element_text(size = 23, 
                                   hjust = 0, 
                                   vjust = 0),
        legend.title = element_text(size = 18),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(5, "mm"),
        axis.title.x = element_text(size = 24),
        axis.title.y =element_text(size = 24,vjust=5),
        plot.title = element_text(size=22,face="bold",hjust=0.5),
        axis.text.x= element_text(size=22,angle=76, vjust=0.73),
        axis.text.y= element_text(size=22),
        strip.background = element_blank(),
        strip.text = element_text(size=23), 
        strip.placement = "outside",
        plot.margin = unit(c(1,1,-2, y.position[fi]), "cm")) 
}

plot.stability.prob.plot <- ggarrange(plotlist=plot.stability.prob,
                           nrow=1,widths=c(1.2,1,1,1),
                           common.legend=T, 
                           legend="bottom", labels = c("A.","B.","C.","D."),
                           label.x = c(0.2,0,0,0),
                           label.y = 0.965,
                           font.label = list(size = 20, 
                                             color = "black", face = "bold",
                                             family = NULL))

plot.stability.prob.plot

ggsave(plot.stability.prob, 
       file = "figures/plot.stability.prob.pdf")

# poster
plot.stability.prob <-  df.stability.summary.small %>%
  mutate(function.name = case_when(function.int==1 ~"1.Traditional",
                                   function.int==2 ~"2.Linear",
                                   function.int==3 ~"3.Exp",
                                   function.int==4 ~"4.Sigmoid")) %>% 
  filter(external_factor =="No external factor") %>%
  filter(comp.com == "two-species community") %>%
  mutate(stab.prob = case_when(stability.significance == "All species" & comp.com =="two-species community"~"2 stable species",
                               stability.significance == "All species" & comp.com =="one-species community"~"1 stable species",
                               stability.significance == "1 species out of 2" ~"1 stable species",
                               stability.significance == "0 species" ~"0 species")) %>%
  aggregate(sim ~ significance + function.name, length) %>%
  mutate(sim = (sim/1500)) %>%
  ggplot() +
  geom_bar(aes(x=function.name, #x=as.factor(stab.prob),
               y=sim, fill=significance ), color="black",
           stat="identity") + 
  labs(x="") +
  #facet_wrap(.~function.name, nrow=1, strip.position= "top") +
  scale_y_continuous("probability of detection in 2 species-communities",
                     labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual("Community dynamics",
                    labels=c("no detected dynamics",
                             "short synchrony","long synchrony",
                             "long & short synchrony",
                             "sychrony & oscilattion"),
                    values=  safe_colorblind_palette_4 ) +
  guides(fill= guide_legend(override.aes = list(pattern = "none"),
                            # position =c(.95, .95),
                            direction="vertical",
                            byrow = TRUE,
                            nrow = 6,
                            title.hjust = 0.1)) +
  theme_minimal() +
  theme(panel.spacing.x = unit(10, "mm"),
        legend.box.background = element_rect(color="black", 
                                             size=0.5),
        legend.position = c(0.34, 0.9),
        panel.grid.major.x  = element_blank(),
        legend.text = element_text(size = 19, 
                                   hjust = 0, 
                                   vjust = 0),
        legend.title = element_text(size = 20),
        legend.key.width = unit(10, "mm"),
        legend.key.height = unit(5, "mm"),
        axis.title = element_text(size = 24),
        axis.text.x= element_text(size=20,angle=76, vjust=0.73),
        axis.text.y= element_text(size=22),
        strip.background = element_blank(),
        strip.text = element_text(size=23), 
        strip.placement = "outside",
        plot.margin = unit(c(1,0,-1,1), "cm")) 
plot.stability.prob

#all
plot.stability.prob.all <-  df.stability.summary.small %>%
  mutate(stab.prob = case_when(stability.significance == "All species" & comp.com =="two-species community"~"two species stable",
                               stability.significance == "All species" & comp.com =="one-species community"~"one species stable",
                               stability.significance == "1 species out of 2" ~"one species stable",
                               stability.significance == "0 species" ~"no species stavle")) %>%
  aggregate(sim ~ stab.prob + significance + function.name + external_factor, length) %>%
  mutate(sim = (sim/1500)) %>%
  ggplot() +
  geom_bar(aes(x=as.factor(stab.prob), y=sim, fill=significance ), color="black",
           stat="identity") + 
  facet_wrap(external_factor~function.name, 
             nrow=3, strip.position= "top") +
  scale_y_continuous("probability of detection\nfor a 1 or 2 species-community",
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete("", 
                   labels=c("no species\nstable","one species\nstable\n","two species\nstable\n"))+
  scale_fill_manual("Community\ndynamics",
                    labels=c("no detected\ndynamics","short\nsynchrony","long\nsynchrony",
                             "long & short\nsynchrony","sychrony &\noscilattion"),
                    values=  safe_colorblind_palette_4 ) +
  guides(fill= guide_legend(override.aes = list(pattern = "none"),
                            direction="vertical",
                            byrow = TRUE,
                            nrow = 5,
                            title.hjust = 0.1)) +
  theme_minimal() +
  theme(panel.spacing.x = unit(10, "mm"),
        panel.grid.major.x  = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 16, 
                                   hjust = 0, 
                                   vjust = 0.5),
        legend.title = element_text(size = 18),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(5, "mm"),
        axis.title = element_text(size = 24),
        axis.text.x= element_text(size=16,angle=66, vjust=0.6),
        axis.text.y= element_text(size=20),
        strip.background = element_blank(),
        strip.text = element_text(size=20), 
        strip.placement = "outside",
        plot.margin = unit(c(1,0,0,1), "cm")) 
plot.stability.prob.all

ggsave(plot.stability.prob.all, 
       file = "figures/plot.stability.prob.all.pdf")




library(ggpubr)
ggarrange(plot.persist.prob,plot.stability.prob,
          ncol=2, common.legend=F,widths =c(1,1.5),
          align =("h"))





#---- heat map -----
df.stability.summary.heatmap <- df.stability.summary.small %>%
  mutate(significance = case_when(significance ==  "synchronous oscillation"|
                                    significance =="oscillation" ~ "oscillatory dynamics" ,
                                  
                                  significance ==  "half-stable synchronous oscillation"|
                                    significance ==  "half-stable oscillation"|
                                    significance == "half-stable synchrony"|
                                    significance == "half-stable"~"half-stable dynamics",
                                  
                                  significance == "stable oscillation"|
                                  significance =="stable synchrony"|
                                    significance =="stable"|
                                    significance =="stable synchronous oscillation" ~ "stable dynamics",
                                  TRUE ~ significance
  )) 

levels(as.factor(df.stability.summary.heatmap$significance))

plot_heatmap_ainit_aslope <- df.stability.summary.heatmap %>%
  gather(a_initial.i.i,a_initial.i.j,a_initial.j.i,a_initial.j.j, 
         key="parameter_a_initial",value="a_initial") %>%
  gather(a_slope.i.i,a_slope.i.j,a_slope.j.i,a_slope.j.j, 
         key="parameter_a_slope",value="a_slope") %>%
  mutate(a_slope = case_when(function.int==1~0,
                             T~a_slope)) %>%
  ggplot(aes(x=a_slope,y=a_initial,
             color=as.factor(significance)))+
  scale_color_manual(values = darken(my_cols, amount = .2)) + 
  scale_fill_manual(values = my_cols) +
  geom_point(alpha=0.5) +
  facet_grid(.~function.int, scales="free") +
  theme_bw()
plot_heatmap_ainit_aslope 

ggsave(plot_heatmap_ainit_aslope ,
       file = "figures/plot_heatmap_ainit_aslope.pdf")


plot_heatmap_c_aslope <- df.stability.summary.heatmap%>%
  gather(c.i.i,c.i.j,c.j.i,c.j.j, 
         key="parameter_c",value="c") %>%
  gather(a_slope.i.i,a_slope.i.j,a_slope.j.i,a_slope.j.j, 
         key="parameter_a_slope",value="a_slope") %>%
  mutate(a_slope = case_when(function.int==1~0,
                             T~a_slope))%>%
  mutate(c = case_when(function.int==1~0,
                       T~c)) %>%
  ggplot(aes(x=a_slope,y=c,
             color=as.factor(significance)))+ 
  geom_point(alpha=0.5) +
  scale_color_manual(values = darken(my_cols, amount = .2)) + 
  scale_fill_manual(values = my_cols) +
  facet_grid(.~function.int, scales="free") +
  theme_bw()
plot_heatmap_c_aslope

plot_heatmap_N_aslope <- df.stability.summary.heatmap %>%
  gather(Nmax.i.i,Nmax.i.j,Nmax.j.i,Nmax.j.j, 
         key="parameter_Nmax",value="Nmax") %>%
  gather(a_slope.i.i,a_slope.i.j,a_slope.j.i,a_slope.j.j, 
         key="parameter_a_slope",value="a_slope") %>%
  mutate(a_slope = case_when(function.int==1~0,
                             T~a_slope))%>%
  mutate(Nmax = case_when(function.int==1~0,
                          T~Nmax)) %>%
  ggplot(aes(x=a_slope,y=Nmax,
             color=as.factor(significance)))+ 
  geom_point(alpha=0.5) +
  scale_color_manual(values = darken(my_cols, amount = .2)) + 
  scale_fill_manual(values = my_cols) +
  facet_grid(.~function.int, scales="free") +
  theme_bw()
plot_heatmap_N_aslope

plot_heatmap_all<- df.stability.summary.heatmap %>%
  gather(Nmax.i.i,Nmax.i.j,Nmax.j.i,Nmax.j.j, 
         key="parameter_Nmax",value="Nmax") %>%
  gather(a_initial.i.i,a_initial.i.j,a_initial.j.i,a_initial.j.j, 
         key="parameter_a_initial",value="a_initial") %>%
  gather(c.i.i,c.i.j,c.j.i,c.j.j, 
         key="parameter_c",value="c") %>%
  gather(a_slope.i.i,a_slope.i.j,a_slope.j.i,a_slope.j.j, 
         key="parameter_a_slope",value="a_slope") %>%
  mutate(a_slope = case_when(function.int==1~0,
                             T~a_slope))%>%
  mutate(Nmax = case_when(function.int==1~0,
                          T~Nmax)) %>%
  mutate(c = case_when(function.int==1~0,
                       T~c)) %>%
  gather(a_initial, a_slope, Nmax, c, 
         key="parameter",value="value") %>%
  ggplot(aes(x=value,y=parameter,
             color=as.factor(significance)))+ 
  geom_point(alpha=0.5) +
  scale_color_manual(values = darken(my_cols, amount = .2)) + 
  scale_fill_manual(values = my_cols) +
  facet_wrap(.~function.int, 
             scales= "free_x") +
  theme_bw()
plot_heatmap_all

ggsave(plot_heatmap_all,
       file = "figures/plot_heatmap_all.pdf")

#---- Sensitivity analysis -----
levels(as.factor(df.stability.summary.heatmap$significance))
df.stability.summary.heatmap$significance <- factor(df.stability.summary.heatmap$significance,
                                                       levels=c("","oscillatory dynamics",
                                                                "half-stable dynamics","stable dynamics"))
df.stability.summary.heatmap$function.int <- factor(df.stability.summary.heatmap$function.int,
                                                    levels=c("1","2",
                                                             "3","4"))
df.stability.summary.sensitivity <- df.stability.summary.heatmap %>%
  mutate(function.int = factor(df.stability.summary.heatmap$function.int,
                               levels=c("1","2","3","4")),
         ratio.i = case_when(is.na(ratio.i) ~ 0,
                             is.infinite(ratio.i) ~ 0,
                             ratio.i > 2 ~ 2,
                             T ~ ratio.i))

model.0 <- lm(ratio.i  ~ function.int,
              df.stability.summary.sensitivity)
model.1 <- lm(formula(paste0("ratio.i  ~ ",
                             paste(c(paste0("a_initial",c(".i.j",".i.i"))),collapse = "+"))), 
              data = df.stability.summary.sensitivity[which(df.stability.summary.sensitivity$function.int==1),])


model.2 <- lm(formula(paste0("ratio.i  ~ ",
                              paste(c(paste0("Nmax",c(".i.j",".i.i")),
                                      paste0("a_initial",c(".i.j",".i.i")),
                                      paste0("a_slope",c(".i.j",".i.i"))),collapse = "+"))), 
               data = df.stability.summary.sensitivity[which(df.stability.summary.sensitivity$function.int==2),])

model.3 <- lm(formula(paste0("ratio.i  ~ ",
                             paste(c(paste0("Nmax",c(".i.j",".i.i")),
                                     paste0("a_initial",c(".i.j",".i.i")),
                                     paste0("c",c(".i.j",".i.i")),
                                     paste0("a_slope",c(".i.j",".i.i"))),collapse = "+"))), 
              data = df.stability.summary.sensitivity[which(df.stability.summary.sensitivity$function.int==3),])

model.4 <- lm(formula(paste0("ratio.i  ~ ",
                             paste(c(paste0("Nmax",c(".i.j",".i.i")),
                                     paste0("a_initial",c(".i.j",".i.i")),
                                     paste0("c",c(".i.j",".i.i")),
                                     paste0("a_slope",c(".i.j",".i.i"))),collapse = "+"))), 
              data = df.stability.summary.sensitivity[which(df.stability.summary.sensitivity$function.int==4),])



sens_out <- tidy(model.4) %>% 
  mutate(function.int = 4) %>% 
  bind_rows(
    (tidy(model.0) %>%
      mutate(term = "function",
             function.int = c(1,2,3,4))),
    (tidy(model.1) %>%
      mutate(function.int = 1)),
    (tidy(model.3) %>%
       mutate(function.int = 3)),
    (tidy(model.2) %>%
       mutate(function.int = 2))
  ) %>%
  filter(term != "(Intercept)")

sens_out$term <- factor(sens_out$term, 
                        levels = 
                          c("function",
                            paste0("a_initial",c(".i.j",".i.i")),
                            paste0("a_slope",c(".i.j",".i.i")),
                            paste0("Nmax",c(".i.j",".i.i")),
                            paste0("c",c(".i.j",".i.i"))
                          )
)
# change to have the focal species first
x_axis_labs <- c("function",
                 paste0("a_initial",c(".i.j",".i.i")),
                 paste0("a_slope",c(".i.j",".i.i")),
                 paste0("Nmax",c(".i.j",".i.i")),
                 paste0("c",c(".i.j",".i.i"))
)


my_cols <- c("#AA4499", "#DDCC77","#88CCEE", "#44AA99")


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
        legend.position = c(.87,.90),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.key.size = unit(0.1, units = "in"),
        axis.text.x = element_text(size = 12, angle = 45, 
                                   hjust = 1, margin = margin(t =10))) +
  #facet_wrap(~comp, nrow = 2) +
  labs(title="Effect size of parameters for each function\n on the stability of species i",
       subtitle= "a_initial.i.i (intra) always > a_initial.i.j (inter)",
       x = "", y = "Effect Size on Stability of i \n ratio mean/var", fill = "function", color = "function") +
  scale_x_discrete(labels = x_axis_labs)

int_sensitivity_plot
ggsave("figures/sensitivity_analysis_stability.pdf", 
       width = 10, height = 8)

#---- detailed visualisation -----
df.detailed.stability <- df.stability.summary %>%
filter(external_factor =="No external factor") %>%
  gather(ratio.i,ratio.j, key="stability",value="value.stability" ) %>%
  filter(!is.na(value.stability)) %>%
  mutate(value.stability =case_when(value.stability > 100 ~ 100,
                                    T ~ value.stability)) %>%
  ggplot(aes(x=comp.com, y=value.stability)) + 
  stat_boxplot() +
  facet_wrap(.~function.name, nrow=1) + 
  theme_bw() + 
  theme(axis.text.x= element_text(size=16,angle=70, hjust=1)) 
  
df.detailed.synchrony <- df.stability.summary %>%
  filter(external_factor =="No external factor") %>%
  filter(!is.na(synchrony.short)) %>%
  ggplot(aes(x=comp.com, y=synchrony.short)) + 
  stat_boxplot() +
  facet_wrap(.~function.name, nrow=1) + 
  theme_bw() + 
  theme(axis.text.x= element_text(size=16,angle=70, hjust=1)) 

df.detailed.oscillation <-  df.stability.summary %>%
    filter(external_factor =="No external factor") %>%
    gather(eigB.j , eigB.i, key= "eigB", value="value.eigB") %>%
    filter(!is.na(value.eigB)) %>%
    ggplot(aes(x=comp.com, y=value.eigB)) + 
    stat_boxplot() +
    facet_wrap(.~function.name, nrow=1) + 
    theme_bw() + 
    theme(axis.text.x= element_text(size=16,angle=70, hjust=1)) 

ggsave(last_plot(),
       file = "figures/summary.abundances.pdf")






