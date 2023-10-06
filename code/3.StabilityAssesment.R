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
 source("code/TimeSerie_toolbox.R")
df.stability <- df.sim.std[which(df.sim.std$invader =="both"),]
 
 summary.df.stability <- NULL
 nsims <- 100
 t.num = 100 # number of generation
 df.stability.summary <- NULL
for(i in 1:nsims){
   for( function.int in 1:4){
     for(add_external_factor in c("none","season","noise")){
       print(paste0("int ", i,"for funct ",function.int, add_external_factor))
       
   df.stability.n.full <-  df.stability[which(df.stability$sim== i &
                                           df.stability$function.int == function.int &
                                           df.stability$external_fact ==add_external_factor),] 
 
   df.stability.n <- df.stability.n.full[c(10:t.num),] # burn first 10 generations
  
#Mean  over time   
   mean.i=c(mean(df.stability.n$Ni))
   mean.j=c(mean(df.stability.n$Nj))
   
 #Mean Variance over time
 # Plaza et al,2012 https://doi.org/10.1111/j.1654-1103.2011.01381.x
 #  amplitude of population fluctuations by means of the standard deviation
   # ratio between median and interquantile
   if(exists("stability.significance")){rm("stability.significance")}
   
   if(!is.na( mean.i) & !is.na( mean.j)){
   ratio_i <- mean(df.stability.n$Ni,na.rm=T)/var(df.stability.n$Ni,na.rm=T) 
   ratio_j <- mean(df.stability.n$Nj,na.rm=T)/var(df.stability.n$Nj,na.rm=T)
   
      msdi <- sd(log10(df.stability.n$Ni))
      msdj <- sd(log10(df.stability.n$Nj))
      

    if(!is.na(ratio_j) & !is.na(ratio_i)){
     if(is.infinite(ratio_i) & is.infinite(ratio_j) ){
        stability.significance = "stable"
         }else{if( ratio_i > 1 & ratio_j > 1){
           stability.significance = "stable"
               }else{if( ratio_i > 1 | ratio_j > 1){
                 stability.significance = "half-stable"
                      }else{stability.significance = "oscillatory"
                            }
                    }
              }
        }
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
      if(mean(log(df.stability.n$Nj))>0 & !is.na(mean(log(df.stability.n$Nj)))){
      fitARIMA.j <- auto.arima(log(df.stability.n$Nj)) # function in R uses a combination of unit root tests, minimization of the AIC and MLE to obtain an ARIMA model
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
        if(mean(log(df.stability.n$Ni))>0 & !is.na(mean(log(df.stability.n$Ni)))){
      fitARIMA.i <- auto.arima(log(df.stability.n$Ni)) 
      p.i = summary(fitARIMA.i)$arma[1]#ARIMA(p,d,q)
      q.i = summary(fitARIMA.i)$arma[2]
      
      if(p.i > 0 & q.i >= 0){
        coeff.i = summary(fitARIMA.i)$coef
        arima.extract.i  <-ARMApqREMLfunct(coeff.i ,log(df.stability.n$Ni),p.i,q.i)
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
     if(aggresShort>1 & aggresLong> 1){
       synchrony.significance <- 1
     }else{synchrony.significance <- 0}
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
                                     msd.i = c(msdi),
                                     msd.j = c(msdj),
                                     ratio.i = c(ratio_i),
                                     ratio.j = c(ratio_j),
                                     oscillation.significance =   oscillation.significance, 
                                     stability.significance = stability.significance,
                                     #correlation.estimate=correlation.estimate,
                                     #correlation.pvalue=correlation.pvalue,
                                     #correlation.significance= correlation.significance,
                                     synchrony.significance = synchrony.significance,
                                     synchrony.long = aggresLong,
                                     synchrony.short = aggresShort) %>%
       bind_cols(bind_cols(coeff.j,coeff.i)) %>%
       bind_cols(df.stability.n.full[1,c(stand.variable,"com.comp.GRWL","com.comp.coex.int")])
    
    df.stability.summary <- bind_rows(df.stability.summary,df.stability.summary.n)
                                     
 
     }
   }
}

save(file="results/df.stability.summary.csv.gz",
     df.stability.summary) 
 
load("results/df.stability.summary.csv.gz") 
#---- global visualisation -----
df.stability.summary.small <- df.stability.summary %>%
  mutate(#correlation.significance = case_when(correlation.significance == 1~"corr",
          #                                    TRUE ~ "no"),
         synchrony.significance = case_when(synchrony.significance == 1~"synchrony",
                                            TRUE ~ "no"))

df.stability.summary.small[df.stability.summary.small =="no"] <- NA

df.stability.summary.small$significance <- apply(df.stability.summary.small[, c("synchrony.significance",
                                              "stability.significance","oscillation.significance")],
                         1, function(x) toString(na.omit(x)))
df.stability.summary.small <- df.stability.summary.small %>%
  mutate(significance = case_when(significance == "synchrony, oscillatory, oscillation" |
                                    significance == "synchrony, oscillatory"|
                                     significance == "synchrony, oscillatory, half_oscillation" ~ "synchronous oscillation",
                                  
                                   significance == "synchrony, half-stable, half_oscillation"|
                                    significance == "synchrony, half-stable, oscillation" ~ "half-stable synchronous oscillation",
                                 
                                   significance =="stable, half_oscillation"|
                                   significance =="stable, oscillation" ~ "stable oscillation",
                                   significance =="oscillatory"|
                                    significance =="oscillatory, half_oscillation" ~ "oscillation",
                                  significance == "half-stable, half_oscillation"|
                                   significance == "half-stable, oscillation" ~ "half-stable oscillation",
                                  significance == "synchrony, half-stable" ~ "half-stable synchrony",
                                  significance == "synchrony, stable" ~ "stable synchrony",
                                  significance == "synchrony, stable, half_oscillation"|
                                    significance == "synchrony, stable, oscillation"~ "stable synchronous oscillation",
                                  TRUE ~ significance
                                  )) 

levels(as.factor(df.stability.summary.small$significance))
levels(as.factor(df.stability.summary.small$external_factor))

  df.stability.summary.small[which(is.na(df.stability.summary.small$significance)),]

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
colorblind_palette_8 <- safe_colorblind_palette[!safe_colorblind_palette %in% my_cols]
#scales::show_col(safe_colorblind_palette)

summary.stability.plot <-  df.stability.summary.small %>%
  filter(com.comp.GRWL == "j_i") %>%
  ggplot(aes(x=as.factor(com.comp.GRWL), fill=as.factor(significance))) + 
  scale_fill_manual(values=  safe_colorblind_palette ) +
  stat_count(position="stack")+
  labs(title ="Percentage of community predicted to have both species in community with underlying dynamics",
      subtitle = " initial intraspecific interactions > initial interspecific community",
      pattern= "Community composition",
      fill = "Community dynamics") + 
  guides(fill= guide_legend(override.aes = list(pattern = c("none")))) +
  facet_wrap(as.factor(external_factor)~function.int, ncol=4, nrow=3) +
  theme(panel.background = element_blank(),
        legend.key.size = unit(1, 'cm')) 
  #theme_bw()

summary.stability.plot 

ggsave(summary.stability.plot, 
       file = "figures/summary.stability.plot.pdf")
#values= c(""="#888888",
#          "extinction"="#661100","half.stable"="#CC6677","stable"="#AA4499",
#          "half_oscillation"="#999933", "oscillation"= "#117733",              
#          "synchrony"="#332288","synchrony, extinction"="#88CCEE",      
#          "synchrony, half_oscillation"="#6699CC", "synchrony, oscillation"="#44AA99"))
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
df.stability.summary <- df.stability.summary
df.stability.summary %>%
  gather(msd.i,msd.i, key="msd",value="value.msd" ) 
  gather(msd.i,msd.i, key="msd",value="value.msd" ) %>%
    df.stability.summary%>%
    ggplot(aes(fill=as.factor(function.int))) + 
           stat_count(aes(x=as.factor(correlation.significance))) +
  theme_bw() + 
  scale_fill_colorblind()
  
  df.stability.summary%>%
    ggplot(aes(fill=as.factor(function.int))) + 
    stat_count(aes(x=synchrony.significance), na.rm= T) +
    theme_bw() + 
    scale_fill_colorblind()

  df.stability.summary%>%
    gather(eigB.j , eigB.i, key= "eigB", value="value.eigB") %>%
    ggplot(aes(x=as.factor(function.int))) + 
    geom_violin(aes(y=value.eigB,color=eigB), na.rm = T) +
    theme_bw() + 
    scale_fill_colorblind()

ggsave(last_plot(),
       file = "figures/summary.abundances.pdf")






