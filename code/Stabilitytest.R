library(forecast) # for ARIMA function
library(lmtest) # for ARIMA eigen value
library(ppcor) # for variable correlation
library(tsvr) # for synchrony 
library(tidyverse) # for synchrony 
library(ggplot2) 
library(ggthemes)
library(ggpattern)
library(wesanderson)
##########################################################################################################
# Stability metric
##########################################################################################################
 source("code/TimeSerie_toolbox.R")
df.stability <- df.sim.std[which(df.sim.std$invader =="both"),]
 
 summary.df.stability <- NULL
 nsims <- 500
 t.num = 100 # number of generation
 df.stability.summary <- NULL
 for(i in 1:nsims){
   for( function.int in 1:4){
     print(paste0("int ", i,"for funct ",function.int))
   df.stability.n <-  df.stability[which(df.stability$sim== i &
                                           df.stability$function.int == function.int),]
 
   df.stability.n <- df.stability.n[c(10:t.num*10),] # burn first 10 generations
  
#Mean  over time   
   mean.i=c(mean(df.stability.n$Ni))
   mean.j=c(mean(df.stability.n$Nj))
   
 #Mean Variance over time
 # Plaza et al,2012 https://doi.org/10.1111/j.1654-1103.2011.01381.x
 #  amplitude of population fluctuations by means of the standard deviation
   # ratio between median and interquantile
   if(exists("stability.significance")){rm("stability.significance")}
   
   if(!is.na( mean.i) & !is.na( mean.j)){
   ratio_i <- median(df.stability.n$Ni,na.rm=T)/IQR(df.stability.n$Ni,na.rm=T) 
   ratio_j <- median(df.stability.n$Nj,na.rm=T)/IQR(df.stability.n$Nj,na.rm=T)
   
      msdi <- sd(log10(df.stability.n$Ni))
      msdj <- sd(log10(df.stability.n$Nj))
      

    if(!is.na(ratio_j) & !is.na(ratio_i)){
     if(is.infinite(ratio_i) & is.infinite(ratio_j) ){
        stability.significance = "stable"
         }else{if( ratio_i > 1 & ratio_j > 1){
           stability.significance = "stable"
               }else{if( ratio_i > 1 | ratio_j > 1){
                 stability.significance = "half-stable"
                      }else{stability.significance = "extinction"
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
        
        
#CORRELATION OF NI AND NJ high value shows that hat x and y are highly consistent and they increase with each other  
      if(exists("correlation.list")){rm("correlation.list")}
      if(!is.na(mean( df.stability.n$Ni)) & !is.na(mean( df.stability.n$Nj))){
      an.error.occured <- FALSE
      tryCatch( {  correlation.list  <-   pcor(data.frame(Ni= df.stability.n$Ni, Nj = df.stability.n$Nj),
                                               method = c("pearson")); print(res) }
                , error = function(e) {an.error.occured <<- TRUE})
      if(!exists("correlation.list")){ 
        correlation.list<- pcor(data.frame(Ni= log(df.stability.n$Ni), Nj = log(df.stability.n$Nj)),
                                method = c("pearson"))
      }
     correlation.estimate <-correlation.list$estimate[1,2] 
     correlation.pvalue <- correlation.list$p.value[1,2]  
    
     if(!is.na(correlation.pvalue) & correlation.pvalue < 0.05 ){
      correlation.significance <- 1
     }else{
       correlation.estimate <- NA
       correlation.pvalue <- NA
       correlation.significance <- NA}
      }else{
        correlation.estimate <- NA
        correlation.pvalue <- NA
        correlation.significance <- NA
      }
     
#SYNCHRONIE OF NI AND NJ - > 1 synchrony | < 1 asynchrony
     if(   msdi > 0.05  & msdj > 0.05 &!is.na(msdi) &!is.na(msdj) &
           mean(df.stability.n$Ni) > 0.05  & mean(df.stability.n$Nj) > 0.05 ){ # otherwise thre are no variation in the data
     df.stability.n.small <- t(as.matrix(data.frame(Ni= df.stability.n$Ni,
                                                    Nj =df.stability.n$Nj)[c(10:nrow( df.stability.n)),]))

     vr.trial <- tsvreq_classic(df.stability.n.small)
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
     
     df.stability.summary.n <- data.frame(sim = as.integer(i),
                                     function.int = as.integer(function.int),
                                     mean.i=c(mean(df.stability.n$Ni)),
                                     mean.j=c(mean(df.stability.n$Nj)),
                                     msd.i = c(msdi),
                                     msd.j = c(msdj),
                                     ratio.i = c(ratio_i),
                                     ratio.j = c(ratio_j),
                                     oscillation.significance =   oscillation.significance, 
                                     stability.significance = stability.significance,
                                     correlation.estimate=correlation.estimate,
                                     correlation.pvalue=correlation.pvalue,
                                     correlation.significance= correlation.significance,
                                     synchrony.significance = synchrony.significance,
                                     synchrony.long = aggresLong,
                                     synchrony.short = aggresShort) %>%
       bind_cols(bind_cols(coeff.j,coeff.i))
    
    df.stability.summary <- bind_rows(df.stability.summary,df.stability.summary.n)
                                     
 
   }
}

save(file="results/df.stability.summary.csv.gz",
     df.stability.summary) 
 
load("results/df.stability.summary.csv.gz") 
#---- global visualisation -----
df.sim.std.small[which(df.sim.std.small$significance==""),]
df.sim.std.small <- df.sim.std[which(df.sim.std$invader == "both" &
                                       df.sim.std$time == 1 #&
                                       #(df.sim.std$function.int == 1| 
                                       #   df.sim.std$function.int == 4)
                                     ),] %>%
  dplyr::select(c("sim","function.int","com.comp.coex", "com.comp.coex.int")) %>%
  left_join(df.stability.summary[,c("sim","function.int","correlation.significance", "synchrony.significance",
                                    "stability.significance","oscillation.significance")], 
            by = c("sim","function.int" )) %>%
  mutate(correlation.significance = case_when(correlation.significance == 1~"corr",
                                              TRUE ~ "no"),
         synchrony.significance = case_when(synchrony.significance == 1~"synchrony",
                                            TRUE ~ "no"))

df.sim.std.small[df.sim.std.small =="no"] <- NA

df.sim.std.small$significance <- apply(df.sim.std.small[, c("synchrony.significance",
                                              "stability.significance","oscillation.significance")],
                         1, function(x) toString(na.omit(x)))
df.sim.std.small <- df.sim.std.small %>%
  mutate(significance = case_when(significance == "synchrony, extinction"|
                                    significance == "synchrony, extinction, half_oscillation"|
                                     significance == "synchrony, extinction, oscillation" ~ "extinction",
                                  #significance == "synchrony, half-stable oscillation"|
                                 #   significance == "synchrony, half-stable oscillation, half_oscillation"|
                                 #   significance == "synchrony, half-stable oscillation, oscillation" ~ "synchrony, half-stable oscillation",
                                #  significance == "synchrony, stable oscillation, half_oscillation"~ "synchrony, stable oscillation",
                                 # significance == "half-stable oscillation, half_oscillation" ~ "half-stable oscillation",
                                 # significance == "stable oscillation, oscillation" ~ "stable oscillation",
                                  TRUE ~ significance
                                  ))

levels(as.factor(df.sim.std.small$significance))
  df.sim.std.small[which(is.na(df.sim.std.small$significance)),]

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
#scales::show_col(safe_colorblind_palette)

summary.stability.plot <-  ggplot(df.sim.std.small, aes(x=as.factor(com.comp.coex),
                                                        fill=as.factor(significance))) + 
  geom_bar_pattern(aes(pattern=as.factor(com.comp.coex.int)),
                   position = "stack",
                   stat= "count",
                   color = "black",
                   pattern_angle =45,
                   pattern_density = 0.2,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 1,
                   na.rm=T) +  
  scale_fill_manual(values=  c("#888888",
                               "#661100","#CC6677","#882255","#AA4499",
                               "#999933",  "#117733",              
                               "#332288","#88CCEE",      
                               "#6699CC", "#44AA99",
                               "#DDCC77")) +
  scale_pattern_manual(values=c("stripe","none"),
                       labels=c("Less than 2 species","2 species")) +
  labs(title ="Percentage of community predicted to have both species in community with underlying dynamics",
      subtitle = " initial intraspecific interactions > initial interspecific community",
      pattern= "Community composition",
      fill = "Community dynamics") + 
  facet_wrap(function.int~., ncol=2, nrow=2) +
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






