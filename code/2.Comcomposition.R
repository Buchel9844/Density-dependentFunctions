library(IDPmisc) # for function NaRV.omit()
##########################################################################################################
# 1. minimum  population  size
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
for(i in 1:nsims){
  for( function.int in 1:4){
    for(add_external_factor in c("No external factor","Noisy change","Periodic change")){
      print(paste0("int ", i,"for funct ",function.int, add_external_factor))
     vec.dNi <- df.sim[which(df.sim$sim.i == i &
                              df.sim$function.int == function.int &
                              df.sim$invader =="Ni" &
                              df.sim$external_factor == add_external_factor),"dNi"]
     vec.dNi.both <- df.sim[which(df.sim$sim.i == i &
                               df.sim$function.int == function.int &
                               df.sim$invader =="both" &
                               df.sim$external_factor == add_external_factor),"dNi"]
     vec.dNj.both <- df.sim[which(df.sim$sim.i == i &
                                    df.sim$function.int == function.int &
                                    df.sim$invader =="both" &
                                    df.sim$external_factor == add_external_factor),"dNj"]
     
     
     vec.dNj <- df.sim[which(df.sim$sim.i == i &
                                df.sim$function.int == function.int &
                                df.sim$invader =="Nj" &
                                df.sim$external_factor == add_external_factor),"dNj"]
     
     
     
     df.min.exp.abun.n <- data.frame(sim = i,
                function.int= function.int,
                external_factor = add_external_factor,
                invader = c("both","both","Ni","Nj"),
                focal = c("Ni","Nj","Ni","Nj"),
                min.exp.abun = c(min.exp.abun(dN=vec.dNi.both,Time = 100, n=0.05),
                      min.exp.abun(dN=vec.dNj.both,Time = 100, n=0.05),
                      min.exp.abun(dN=vec.dNi,Time = 100, n=0.05),
                      min.exp.abun(dN=vec.dNj,Time = 100, n=0.05)))
     df.min.exp.abun <- bind_rows(df.min.exp.abun,df.min.exp.abun.n)
    }
  }
}
 
df.min.exp.abun <- df.min.exp.abun %>%
  mutate(function.name = case_when(function.int==1 ~"1.Constant",
                                   function.int==2 ~"2.Linear",
                                   function.int==3 ~"3.Exp",
                                   function.int==4 ~"4.Sigmoid"))
         
write.csv(df.sim , 
          file = paste0("results/df.sim_add_external_factor.csv.gz"))

##########################################################################################################
# 1. Compute slopes of growth over time
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



