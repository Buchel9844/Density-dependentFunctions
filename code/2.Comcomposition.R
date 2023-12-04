library(IDPmisc) # for function NaRV.omit()
library(ggpattern)
##########################################################################################################
# 1. minimum  abundance
##########################################################################################################
df.sim <- as.data.frame(df.sim)
df.min.abundance <- NULL
nsims <- 750
for(i in 1:nsims){
  for( function.int in 1:4){
    for(add_external_factor in c("No external factor")){
      print(paste0("int ", i,"for funct ",function.int, add_external_factor))
      df.min.abundance.n <- df.sim[which(df.sim$sim.i == i &
                        df.sim$function.int == function.int &
                        df.sim$external_factor == add_external_factor &
                          df.sim$time > 99),]

      vec.abundance.Ni <- df.min.abundance.n$Ni
      vec.abundance.Nj <- df.min.abundance.n$Nj
      vec.GR.Ni <- df.min.abundance.n$dNi
      vec.GR.Nj <- df.min.abundance.n$dNj
      
      df.runaway.n <- df.sim[which(df.sim$sim.i == i &
                                           df.sim$function.int == function.int &
                                           df.sim$external_factor == add_external_factor),]
       
  
      df.min.abun.n <- data.frame(sim = i,
                                  a_initial.inter = c(df.min.abundance.n$a_initial.i.j[1],df.min.abundance.n$a_initial.j.i[1]),
                                 function.int= function.int,
                                external_factor = add_external_factor,
                                focal = c("Ni","Nj"),
                                mean.abundance = c(mean(vec.abundance.Ni, na.rm=T),mean(vec.abundance.Nj, na.rm=T)),
                                median.abundance = c(median(vec.abundance.Ni, na.rm=T),median(vec.abundance.Nj, na.rm=T)),
                                min.abundance = c(min(vec.abundance.Ni, na.rm=T),min(vec.abundance.Nj, na.rm=T)),
                                min.GR = c(vec.GR.Ni[which(vec.abundance.Ni ==min(vec.abundance.Ni, na.rm=T))[1]],
                                           vec.GR.Nj[which(vec.abundance.Nj ==min(vec.abundance.Nj, na.rm=T))[1]]),
                                max.abundance = c(max(vec.abundance.Ni, na.rm=T),max(vec.abundance.Nj, na.rm=T)),
                                max.GR = c(max(vec.GR.Ni, na.rm=T),max(vec.GR.Nj, na.rm=T)),
                                max.abundance.all = c(max(df.runaway.n$Ni, na.rm=T),max(df.runaway.n$Nj, na.rm=T)),
                                max.GR.all = c(max(df.runaway.n$dNi, na.rm=T),max(df.runaway.n$dNj, na.rm=T)),
                                median.abundance.all = c(median(df.runaway.n$Ni, na.rm=T),median(df.runaway.n$Nj, na.rm=T)))
      
      df.min.abundance <- bind_rows(df.min.abundance,df.min.abun.n)
    }
  }
}

df.min.abundance  <- df.min.abundance  %>%
  mutate(function.name = case_when(function.int==1 ~"1.Constant",
                                   function.int==2 ~"2.Linear",
                                   function.int==3 ~"3.Exp",
                                   function.int==4 ~"4.Sigmoid")) %>%
  mutate(Mean.class = case_when(mean.abundance <=0 ~"Null",
                                mean.abundance > 500 ~"Runaway",
                              T ~"Okay"),
         Median.class = case_when(median.abundance <=0 ~"Null",
                                  median.abundance.all > 1000 ~"Runaway",
                                T ~"Okay"),
         GR.class = case_when(min.GR >= 1  & max.GR < 500 ~"Positive",
                              max.GR.all > 500 ~"unbound" ,
                              min.GR==0 & max.GR < 500 ~"Null constant",
                                   T ~"Negative or NA"),
         Abundance.class = case_when( min.abundance > 0 & max.abundance < 1000 ~ "Persisting",
                                      max.abundance.all  > 1000 ~"unbound",
                                     T ~"Extinct or NA"))
         
#save dataset
write.csv(df.min.abundance , 
          file = paste0("results/df.min.abundance.csv.gz"))
df.min.abundance = data.table::fread("results/df.min.abundance.csv.gz")
head(df.min.abundance)
str(df.min.abundance)
df.min.abundance <- df.min.abundance[,-1]

# make it horyzontal
df.min.abun.horyzontal <- full_join(df.min.abundance[which(df.min.abundance$focal =="Ni"),],
                                    df.min.abundance[which(df.min.abundance$focal =="Nj"),], suffix = c(".Ni",".Nj"),
          by=c("function.int","sim","external_factor","function.name")) %>%
          mutate(comp.com = case_when((Median.class.Ni =="Null" & Median.class.Nj =="Null" )|
                                       (GR.class.Ni=="Negative or NA" & GR.class.Nj=="Negative or NA") ~"no species",
                                      
                                      GR.class.Ni=="Positive" & GR.class.Nj=="Positive" &
                                      Median.class.Ni =="Okay" & Median.class.Nj =="Okay" ~"two-species community",
                                      
                                      (GR.class.Ni=="unbound" | GR.class.Nj=="unbound"|
                                         Median.class.Ni =="Runaway"|Median.class.Nj =="Runaway"|
                                         Abundance.class.Ni =="unbound"|Abundance.class.Nj =="unbound") ~ "run away population",
                                      
                               T ~"one-species community")) %>%
  mutate(Interspecific.interaction = case_when(a_initial.inter.Ni > 0 &
                                                 a_initial.inter.Nj  > 0 ~ "facilitation",
                                               a_initial.inter.Ni < 0 &
                                                 a_initial.inter.Nj  < 0 ~ "competition",
                                               T~"both")) %>%
  mutate(comp.com = factor(comp.com, levels=c("run away population","no species","one-species community","two-species community")),
         Interspecific.interaction = factor(Interspecific.interaction,
                                            levels=c("competition","both","facilitation")))


# visualisation
#my_cols <- c("#AA4499", "#DDCC77","#88CCEE", "#44AA99")
cols_interaction <- c("#661100", "#888888", "#6699CC", "#332288") # "#DDCC77","#661100","#117733")
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")



com.comp.plot <- ggplot(df.min.abun.horyzontal[df.min.abun.horyzontal$external_factor =="No external factor",],
                            aes( fill=as.factor(comp.com),
                                 x=as.factor(Interspecific.interaction) )) +

  geom_bar()+
  #labels=c("run away population","no species","1-species","2-species")) +
  scale_color_manual(values = darken(cols_interaction, amount = .1)) + 
  #scale_pattern_fill_manual(values = my_cols) + 
  scale_fill_manual(
  values = cols_interaction) + 
  facet_wrap(.~function.name,nrow=1) +
  scale_y_continuous(limits=c(0,250), expand = c(0, 0)) +
 # scale_x_discrete(labels = c("run away","no species","one-species","two-species")) +
  theme_minimal() +
  labs(y="Number of simulated communities", fill="Community trajectory",
       x="Governing interaction")+
       #title="Number of communities with one or two species \nhaving a positive or null growth rate \nwhen low AND a positive abundance")+
  guides(color="none") +
  theme( legend.key.size = unit(1, 'cm'),
         legend.position = "bottom",
         strip.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(),
         strip.text = element_text(size=28),
         legend.text=element_text(size=20),
         legend.title=element_text(size=20),
         #axis.ticks.x=element_blank(),
         axis.text.x= element_text(size=20, angle=66, hjust=1),
         axis.text.y= element_text(size=20),
         axis.title.x= element_text(size=24),
         axis.title.y= element_text(size=24),
         title=element_text(size=16))

com.comp.plot
ggsave(com.comp.plot,
       file = "figures/com.comp.plot.pdf")

# ALL
com.comp.plot.all <- ggplot(df.min.abun.horyzontal,
                            aes( fill=as.factor(comp.com),
                                 x=as.factor(Interspecific.interaction) )) +
  
  geom_bar()+
  #labels=c("run away population","no species","1-species","2-species")) +
  scale_color_manual(values = darken(cols_interaction, amount = .1)) + 
  #scale_pattern_fill_manual(values = my_cols) + 
  scale_fill_manual(
    values = cols_interaction) + 
  facet_wrap(external_factor~function.name) +
  scale_y_continuous(limits=c(0,250), expand = c(0, 0)) +
  # scale_x_discrete(labels = c("run away","no species","one-species","two-species")) +
  theme_minimal() +
  labs(y="Number of simulated communities", fill="Community trajectory",
       x=expression(paste("Initial interaction value", alpha['0,i,j'])))+
  #title="Number of communities with one or two species \nhaving a positive or null growth rate \nwhen low AND a positive abundance")+
  guides(color="none") +
  theme( legend.key.size = unit(1, 'cm'),
         legend.position = "right",
         strip.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(),
         strip.text = element_text(size=20),
         legend.text=element_text(size=16),
         legend.title=element_text(size=16),
         #axis.ticks.x=element_blank(),
         axis.text.x= element_text(size=16, angle=66, hjust=1),
         axis.text.y= element_text(size=20),
         title=element_text(size=16))

com.comp.plot.all
ggsave(com.comp.plot.all,
       file = "figures/com.comp.plot.all.pdf")


# join dataset with full simulation
df.sim$sim <- df.sim$sim.i
df.sim.comcomp <- left_join(df.sim,
                            df.min.abun.horyzontal,
                            by=c("function.int","sim","external_factor","function.name"))


##########################################################################################################
# 2. minimum  population  size
##########################################################################################################
min.exp.abun <- function(dN,Time,n){ 
  # Time is the time given to reach a particular threshold n
  # n is thethreshold  population  size  as  a  proportion  of  the  initial population  size,  so  that  0  ≤n ≤1
  log.vec.dN <- NaRV.omit(log(dN))
  varN <- sqrt(var(log.vec.dN))
  meanN <- mean(log.vec.dN)
  u = meanN*sqrt(Time)/varN
  v = -log(n)/(varN*sqrt(Time))
  Q =  pnorm(c(-u-v), mean = 0, sd = 1) + exp(-2*u*v)* pnorm(c(u-v),mean = 0, sd = 1)
  return(Q)
}
df.min.exp.abun <- NULL
nsims <- 500
for(i in 1:nsims){
  for( function.int in 1:4){
    for(add_external_factor in c("No external factor","Noisy change","Periodic change")){
      print(paste0("int ", i,"for funct ",function.int, add_external_factor))
     
     vec.dNi.both <- df.sim[which(df.sim$sim.i == i &
                               df.sim$function.int == function.int &
                               df.sim$invader =="both" &
                               df.sim$external_factor == add_external_factor),"dNi"]
     vec.dNj.both <- df.sim[which(df.sim$sim.i == i &
                                    df.sim$function.int == function.int &
                                    df.sim$invader =="both" &
                                    df.sim$external_factor == add_external_factor),"dNj"]
     
     vec.dNi <- df.sim[which(df.sim$sim.i == i &
                               df.sim$function.int == function.int &
                               df.sim$invader =="Ni" &
                               df.sim$external_factor == add_external_factor),"dNi"]
     
     vec.dNj <- df.sim[which(df.sim$sim.i == i &
                                df.sim$function.int == function.int &
                                df.sim$invader =="Nj" &
                                df.sim$external_factor == add_external_factor),"dNj"]
     
     
     
     df.min.exp.abun.n <- data.frame(sim = i,
                function.int= function.int,
                external_factor = add_external_factor,
                invader = c("both","both","Ni","Nj"),
                focal = c("Ni","Nj","Ni","Nj"),
                min.exp.abun = c(min.exp.abun(dN=vec.dNi.both,Time = 100, n=0.001),
                      min.exp.abun(dN=vec.dNj.both,Time = 100, n=0.001),
                      min.exp.abun(dN=vec.dNi,Time = 100, n=0.001),
                      min.exp.abun(dN=vec.dNj,Time = 100, n=0.001)))
     df.min.exp.abun <- bind_rows(df.min.exp.abun,df.min.exp.abun.n)
    }
  }
}
 
df.min.exp.abun <- df.min.exp.abun %>%
  mutate(function.name = case_when(function.int==1 ~"1.Constant",
                                   function.int==2 ~"2.Linear",
                                   function.int==3 ~"3.Exp",
                                   function.int==4 ~"4.Sigmoid"))
    
ggplot(df.min.exp.abun,
       aes( y=min.exp.abun, x=invader)) +
  #geom_density() +
  facet_wrap(external_factor~function.name) +
  theme_bw() +
  geom_boxplot()

write.csv(df.min.exp.abun , 
          file = paste0("results/df.min.exp.abun.csv.gz"))

##########################################################################################################
# 3. Compute slopes of growth over time
##########################################################################################################

load("results/df.sim.csv.gz")

df.glm_all <- NULL

for(i in 1:nsims){
  for( function.int in 1:4){
    for(add_external_factor in c("none","season","noise")){
      print(paste0("int ", i,"for funct ",function.int, add_external_factor))
      

    # Slope of Ni with init.cond= c(1,Nj*)
    df.Ni <- df.sim[which(df.sim$sim.i == i &
                                  df.sim$function.int == function.int &
                                  df.sim$invader =="Ni" &
                                  df.sim$external_factor == add_external_factor),]
    if(sum(is.na(df.Ni$Ni))>1 |df.Ni$Nj[1] == 0|df.Ni$Ni[2] >1000 | df.Ni$Nj[2] >1000){
      df.Ni.glm <-  data.frame(com.comp.prob.Ni =F, 
                               sim = i,
                               function.int = function.int,
                               external_factor = add_external_factor,
                               intercept.Ni =NA,
                               time.Ni = NA,
                               com.comp.Ni = 0)
    }else{
    
    df.Ni.glm <- as.data.frame(summary(glm(formula = Ni ~ time, data = df.Ni,
                                           family = "gaussian"))$coefficients) %>%
      rownames_to_column(var="coeff") %>%
      mutate(coeff = case_when(coeff =="(Intercept)" ~ "intercept",
                               coeff =="time"~"time")) %>%
      dplyr::select(coeff,Estimate ) %>%
      mutate(coeff = paste0( coeff,".Ni")) %>%
      spread(coeff,Estimate) %>%
      mutate(com.comp.Ni = mean(intercept.Ni + 100*time.Ni),
             com.comp.prob.Ni = case_when(com.comp.Ni>0~T,
                                          com.comp.Ni<0~F),
             sim = i,
             external_factor = add_external_factor,
             function.int = function.int)
    }
    # Slope of Nj with init.cond= c(Ni*,1)
    df.Nj <- df.sim[which(df.sim$sim.i == i &
                                  df.sim$function.int == function.int &
                                  df.sim$invader =="Nj" &
                            df.sim$external_factor == add_external_factor),]
    
    if(sum(is.na(df.Nj$Nj))>1 |df.Nj$Ni[1] == 0|df.Nj$Nj[2] >1000 | df.Nj$Ni[2] >1000){
      df.Nj.glm <-  data.frame(com.comp.prob.Nj =F, 
                               sim = i,
                               external_factor = add_external_factor,
                               function.int = function.int,
                               intercept.Nj =NA,
                               time.Nj = NA,
                               com.comp.Nj = 0)
    }else{
    df.Nj.glm <- as.data.frame(summary(glm(formula = Nj ~ time, data = df.Nj,
                                           family = "gaussian"))$coefficients) %>%
      rownames_to_column(var="coeff") %>%
      mutate(coeff = case_when(coeff =="(Intercept)" ~ "intercept",
                               coeff =="time"~"time")) %>%
      dplyr::select(coeff,Estimate ) %>%
      mutate(coeff = paste0( coeff,".Nj"))%>%
      spread(coeff,Estimate) %>%
      mutate(com.comp.Nj = mean(intercept.Nj + 100*time.Nj),
             com.comp.prob.Nj = case_when(com.comp.Nj>0~T,
                                          com.comp.Nj<0~F),
             sim = i,
             external_factor = add_external_factor,
             function.int = function.int)
    }
    # Slope of Nj with init.cond= c(1,1)
    df.Ni.both <- df.sim[which(df.sim$sim.i == i &
                            df.sim$function.int == function.int &
                            df.sim$invader =="both" &
                              df.sim$external_factor ==add_external_factor),]
    
    if(sum(is.na(df.Ni.both$Ni))>1 |df.Ni.both$Nj[1] == 0|df.Ni.both$Ni[2] >1000 | df.Ni.both$Nj[2] >1000){
      df.Ni.both.glm <-  data.frame(com.comp.prob.Ni.both =F, 
                               sim = i,
                               external_factor = add_external_factor,
                               function.int = function.int,
                               intercept.Ni.both =NA,
                               time.Ni.both = NA,
                               com.comp.Ni.both = 0)
    }else{
    df.Ni.both.glm <- as.data.frame(summary(glm(formula = Ni ~ time, data = df.Ni,
                                           family = "gaussian"))$coefficients) %>%
      rownames_to_column(var="coeff") %>%
      mutate(coeff = case_when(coeff =="(Intercept)" ~ "intercept",
                               coeff =="time"~"time")) %>%
      dplyr::select(coeff,Estimate ) %>%
      mutate(coeff = paste0( coeff,".Ni.both"))%>%
      spread(coeff,Estimate) %>%
      mutate(com.comp.Ni.both = mean(intercept.Ni.both + 100*time.Ni.both),
             com.comp.prob.Ni.both = case_when(com.comp.Ni.both >0 ~ T,
                                               com.comp.Ni.both < 0 ~ F),
             sim = i,
             external_factor = add_external_factor,
             function.int = function.int)
    }
    # Slope of Nj with init.cond= c(1,1)
  
    df.Nj.both <- df.sim[which(df.sim$sim.i == i &
                            df.sim$function.int == function.int &
                            df.sim$invader =="both" &
                              df.sim$external_factor == add_external_factor),]
    if(sum(is.na(df.Nj.both$Nj))>1 |df.Nj.both$Ni[1] == 0|df.Nj.both$Nj[2] >1000 | df.Nj.both$Ni[2] >1000){
      df.Nj.both.glm <-  data.frame(com.comp.prob.Nj.both =F, 
                               sim = i,
                               external_factor = add_external_factor,
                               function.int = function.int,
                               intercept.Nj.both =NA,
                               time.Nj.both = NA,
                               com.comp.Nj.both = 0)
    }else{
    df.Nj.both.glm <- as.data.frame(summary(glm(formula = Nj ~ time, data = df.Nj,
                                           family = "gaussian"))$coefficients) %>%
      rownames_to_column(var="coeff") %>%
      mutate(coeff = case_when(coeff =="(Intercept)" ~ "intercept",
                               coeff =="time"~"time")) %>%
      dplyr::select(coeff,Estimate ) %>%
      mutate(coeff = paste0( coeff,".Nj.both"))%>%
      spread(coeff,Estimate) %>%
      mutate(com.comp.Nj.both = mean(intercept.Nj.both + 100*time.Nj.both),
             com.comp.prob.Nj.both = case_when(com.comp.Nj.both>0~T,
                                          com.comp.Nj.both<0~F),
             sim = i,
             external_factor = add_external_factor,
             function.int = function.int)
    
    }
    df.glm <- full_join(full_join(df.Ni.glm,df.Nj.glm),
                         full_join(df.Ni.both.glm,df.Nj.both.glm)) %>%
      mutate(com.comp.coex.int = case_when(mean(com.comp.Ni,com.comp.Ni.both) > 0 &
                                         mean(com.comp.Nj,com.comp.Nj.both) > 0 ~ "j_i",
                                       mean(com.comp.Ni,com.comp.Ni.both) <= 0 & 
                                         mean(com.comp.Nj,com.comp.Nj.both) > 0 ~ "j",
                                       mean(com.comp.Ni,com.comp.Ni.both) > 0 & 
                                         mean(com.comp.Nj,com.comp.Nj.both) <= 0 ~ "i",
                                       mean(com.comp.Ni,com.comp.Ni.both) <= 0 & 
                                         mean(com.comp.Nj,com.comp.Nj.both) <= 0 ~ "0"))
    
    df.glm_all <- bind_rows(df.glm_all,df.glm)
    }
  }
}

df.glm_all[which(is.na(df.glm_all$com.comp.SLOPE)),"com.comp.coex"] <- 0
head(df.glm_all)
str(df.glm_all)

save(df.glm_all,
     file="results/df.glm_all.csv.gz")

load("results/df.glm_all.csv.gz")


plot_com.comp.SLOPE <- ggplot(df.glm_all,aes(x=com.comp.coex.int,
                      group=as.factor(function.int),
                      fill=as.factor(function.int)))+ 
  stat_count()+ theme_bw()  +
  scale_fill_colorblind() + 
  facet_wrap(external_factor~.) +
  labs(title ="Community composition based on mean slope")
plot_com.comp.SLOPE
 
ggsave(plot_com.comp.SLOPE,
       file = "figures/com.comp.SLOPE.pdf")

##########################################################################################################
# 2. Compute minimal growth rate
##########################################################################################################

load("results/df.sim.csv.gz")

df.GRWL <- NULL

for(i in 1:nsims){
  for( function.int in 1:4){
    for(add_external_factor in c("none","season","noise")){
      print(paste0("int ", i,"for funct ",function.int, add_external_factor))
    
    # Slope of Ni with init.cond= c(1,Nj*)
    df.GRWL_n <- df.sim[which(df.sim$sim.i == i &
                            df.sim$function.int == function.int &
                              df.sim$external_factor ==add_external_factor),] %>%
      select(-invader) %>%
      mutate_if(is.character,as.numeric) %>%
      mutate(external_factor = add_external_factor )

    
    GRWL_i <- df.GRWL_n$dNi[which(df.GRWL_n$Ni ==min(df.GRWL_n$Ni))]
    GRWL_i <- min(GRWL_i[which(GRWL_i>=1)])

    if(is.na(GRWL_i)){
      GRWL_i = 0
    }
    
    #lower.quantile_j = quantile(df.GRWL_n$Nj, 0.01, na.rm=T)
    GRWL_j <- df.GRWL_n$dNj[which(df.GRWL_n$Nj ==min(df.GRWL_n$Nj))]
    GRWL_j <- min(GRWL_j[which(GRWL_j>=0)])
    
    if(is.na(GRWL_j)){
      GRWL_j = 0
    }
    
    df.GRWL_n <- df.GRWL_n[1,] %>%
      mutate(GRWL_i=GRWL_i,
             GRWL_j=GRWL_j,
             com.comp.GRWL = case_when(GRWL_i >= 1 & GRWL_j >= 1 ~"j_i",
                                       GRWL_i > 1 & GRWL_j < 1 ~"i",
                                       GRWL_i < 1 & GRWL_j > 1 ~"j",
                                       GRWL_i < 1 & GRWL_j < 1 ~"0")
               
             )
    df.GRWL <- bind_rows(df.GRWL,df.GRWL_n)
    }
  }
}
  
plot_com.comp.GRWL <- ggplot(df.GRWL,aes(x=com.comp.GRWL,
                                        group=as.factor(function.int),
                                        fill=as.factor(function.int)))+ 
  stat_count()+ theme_bw()  +
  scale_fill_colorblind() +
  facet_wrap(external_factor~.) +
  labs(title ="Community composition based on GRWlow")
plot_com.comp.GRWL

ggsave(plot_com.comp.GRWL,
       file = "figures/com.comp.GRWL.pdf")

#join both methods of community composition
df.GRWL$sim <- df.GRWL$sim.i
plot_com.comp <- full_join(df.GRWL,df.glm_all) %>%
  gather(com.comp.GRWL,com.comp.coex.int,
         key="com.comp.int",value="com.comp.id") %>%
  filter(com.comp.id=="j_i") %>%
  mutate(com.comp.int = case_when(com.comp.int=="com.comp.GRWL"~ "Growth rate when low",
                                  com.comp.int=="com.comp.coex.int" ~"Mean decay over time")) %>%
  ggplot(aes(x=com.comp.int, group=as.factor(function.int),
                    fill=as.factor(function.int)))+ 
  stat_count(position="dodge") +
  scale_fill_manual(values=my_cols) +
  labs(fill="interaction \nfunctions",
       x="community composition",
       y="number of simulated communities") +
  facet_grid(external_factor~.) +
  theme_bw()

ggsave(plot_com.comp,
       file = "figures/com.comp.pdf")

##########################################################################################################
# 3. Bind all data
##########################################################################################################



