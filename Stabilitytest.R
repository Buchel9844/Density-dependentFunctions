library(forecast) # for ARIMA function
library(lmtest) # for ARIMA eigen value
library(ppcor) # for variable correlation
library(tsvr) # for synchrony 
library(tidyverse) # for synchrony 
library(ggplot2) 
library(ggthemes)
library(ggpattern)
##########################################################################################################
# Stability metric
##########################################################################################################
 df.stability <- df.sim[which(df.sim$invader =="both"),]
 
 summary.df.stability <- NULL
 nsims <- 1000
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
      msdi <- sd(log10(df.stability.n$Ni))
      msdj <- sd(log10(df.stability.n$Nj))
      if(exists("correlation.list")){rm("stability.significance")}
      
      if((msdi==0 | msdj ==0 | is.na(msdi)| is.na(msdj) ) &
         (mean.i==0|mean.j==0| is.na(mean.i)| is.na(mean.j) )){
        stability.significance = "extinction"
      }else{
         if((msdi==0 | msdj ==0 | is.na(msdi)| is.na(msdj)) & 
        ( mean.i!=0|mean.j!=0| !is.na(mean.i)| !is.na(mean.j) )){
        stability.significance = "half.stable"
         }
      }
     if(msdi==0 & msdj ==0 & !is.na(msdi) & !is.na(msdj) &
         mean.i!=0 & mean.j!=0){
        stability.significance = "stable"
      }
      if(!exists("stability.significance")){
        stability.significance = "no"
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
 
stand.variable <- c("mean","median","var","upperbound","lowerbound")


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
    gather(q.j , q.i, key= "q", value="value.q") %>%
    mutate(value.q.simple = case_when(value.q >=1 ~ 1,
                                      TRUE ~ value.q)) %>%
    ggplot(aes(fill=as.factor(function.int))) + 
    stat_count(aes(x=value.q.simple), na.rm = T) +
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

df.sim.std.small <-df.sim.std[which(df.sim.std$invader == "both" &
                                    df.sim.std$time == 1 ),] %>%
  select("sim","function.int","com.comp.coex", "com.comp.coex.int") %>%
  left_join(df.stability.summary[,c("sim","function.int","correlation.significance", "synchrony.significance")], 
            by = c("sim","function.int" )) %>%
  mutate(correlation.significance = case_when(correlation.significance == 1~"corr",
                                              TRUE ~ "no"),
         synchrony.significance = case_when(synchrony.significance == 1~"synchrony",
                                              TRUE ~ "no")) %>%
  gather(correlation.significance, synchrony.significance, key= "state", value="significance")
str(df.sim.std.small)
summary.stability.plot <-  ggplot(df.sim.std.small, aes(x=as.factor(function.int),
                             fill=as.factor(com.comp.coex.int), 
                             pattern=significance)) + 
  geom_bar_pattern(position = "fill",stat= "count",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6,
                   na.rm=T) +
  scale_pattern_manual(values = c("no" = "none", "synchrony"= "stripe","corr" = 'wave'),
                       labels = c("Correlated abundances", "other","Synchrony of abundances")) +  
  scale_fill_manual(values=c("#E69F00", "#56B4E9"),
                    labels=c("Less than 2 species","2 species")) +
  labs(title ="Percentage of community predicted to have both species in community with underlying dynamics",
       subtitle = " initial intraspecific interactions > initial interspecific community",
       fill= "Community composition",
       pattern = "Community dynamics over time") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")))+ 
  theme_bw()

ggsave(summary.stability.plot, 
       file = "figures/summary.stability.plot.pdf")


##########################################################################################################
# Time series analysis
##########################################################################################################


#Analysis of ecological time series with ARMA(p,q) models
#Anthony R. Ives, Karen C. Abbott, Nicolas L. Ziebarth
#First published: 01 March 2010 https://doi.org/10.1890/09-0442.1

ARMApqREMLfunct <- function(par,nX,p,q){  
  
  ##generates LL for ARMApq model.
  ##translated from matlab function ARMApqREMLfunct.m by Nic Ziebarth and Tony Ives
  ##10 March, 2009
  
  ##par should have values for both AR and MA components.
  ##nX is the data, p = # of AR, q = # of MA
  
  ##reoccuring issue is what should diagExtend output when size = 0. This comes up
  ##when p=1 and q=1. I think any 1 x 1 matrix will work. 
  
  ## this is standardized so that a(1)=1 and is not estimated
  b <- par[1:p];
  lengthPar <- length(par);
  
  ##this slight change relative to Matlab case. If q=0 in old version, then a =[], which ##was no good
  if(q > 0){ 
    a <- c(1,par[(p+1):lengthPar]);
  }
  else{
    a <- 1
  }
  
  B <- diagExtend(p,-1);
  sizeB <- dim(B);
  B[1,1:p] <- b;
  eigs <- eigen(B)$values;
  eigB <- max(abs(eigs));
  
  if (eigB>=1 | max(abs(a))>10){
    LL<-10^10;
  }
  
  T <- length(nX);
  
  ## solve for stationary distribution
  
  k <- max(p,q);
  
  ##deal with case of no MA terms separately
  if (q > 0){
    A <- array(0,dim=c(k,q));
    A[1,1:q] <- a[2:(q+1)];
    if(k > 1){
      B <- diagExtend(k,-1);
      B[1,1:p] <- b;
    }
    else{
      B <- array(1,dim=c(1,1))*b;
    }
    
    CC <- cbind(B,A);
    
    otherPart <- rbind(array(0,dim=c(k,q)), diagExtend(q,1));
    CC <- rbind(CC,t(otherPart));
    
    Ve <- array(0,dim=c(k+q,k+q));
    Ve[1,1] <- 1;
    Ve[1,k+1] <- 1;
    Ve[k+1,1] <- 1;
    Ve[k+1,k+1] <- 1;
    
    sizeVe <- dim(Ve);
    vecVe <- array(Ve,dim=c(sizeVe[1]*sizeVe[2],1));
    C <- solve(diag(1,(k+q)^2) - kronecker(CC,CC))%*%vecVe;
    
    C <- array(C,dim=c(k+q,k+q));
    
    Vstationary<-C[1:k,1:k];
  }
  else{
    
    ##this is particular example of issue when p=1
    if(p > 1){
      B <- diagExtend(p,-1);
      sizeB <- dim(B);
      B[1,1:sizeB[2]] <- b;
    }
    else{
      B <- array(1,dim=c(1,1))*b;
    }
    
    Ve <- array(0,dim=c(p,p));
    Ve[1,1] <- 1;
    
    sizeVe <- dim(Ve);
    vecVe <- array(Ve,dim=c(sizeVe[1]*sizeVe[2],1));	
    C <- solve(diag(1,(p+q)^2) - kronecker(B,B))%*%vecVe;
    
    Vstationary <- array(C,dim=c(p+q,p+q));
  }
  
  ## solve for V_T
  
  ##MA component
  
  ##This was added to separately deal with case of no MA lags
  if(q > 0){
    AA <- diag(1,q+1);
    for (i in 1:q){
      AA <- AA + a[i+1]*diagExtend(q+1,-i);
    }
    AA <- t(AA)%*%AA;
    sizeAA <- dim(AA);
    
    aa <- AA[1,1:sizeAA[2]];
    AA <- array(0,dim=c(T,T));
    for (i in 1:q){
      AA <- AA + aa[i+1]*diagExtend(T,-i);
    }
    
    AA <- AA+t(AA)+aa[1]*diag(1,T);
  }
  else{
    AA <- diag(1,T);
  }
  
  BB <- diag(1,k);
  if(k > 1){ 
    for (i in 1:p){
      BB <- BB - b[i]*diagExtend(k,-i);
    }
  }
  A1 <- BB%*%Vstationary%*%t(BB);
  AA[1:k,1:k] <- A1;
  
  ## AR component
  W <- diag(1,T);
  for (i in 1:p){
    W <- W - b[i]*diagExtend(T,-i);
  }
  invW <- solve(W);
  
  
  V <- invW%*%AA%*%t(invW);
  
  ## compute LL function for extant data
  
  pick <- which(is.na(nX)==FALSE);
  T <- length(pick);
  X <- nX[pick];
  V <- V[pick,pick];
  
  invV <- solve(V);
  
  U <- array(1,dim=c(T,1));
  mu <- solve(t(U)%*%invV%*%U)%*%t(U)%*%invV%*%X;
  H <- X - mu; 
  
  ## condensed ML LL
  ##s2 <- (t(H)%*%invV%*%H)/T;
  ##LL <- .5*(T*log(2*pi)+T*log(s2)+log(det(V))+T);
  
  ## condensed REML LL
  s2 <- (t(H)%*%invV%*%H)/(T-1);
  LL <- .5*((T-1)*log(2*pi)+(T-1)*log(s2)+log(det(V))+log(det(t(U)%*%invV%*%U))+(T-1));
  
  
  if(abs(Im(LL)) > 10^-6) {
    LL <- 10^10;
  }
  return(list(LL=LL,eigB=eigB))
}

## end ARMApqREMLfunct
diagExtend <- function(p,diagInd){
  
  ##This is a function with the functionality of diag in Matlab
  ##In particular, this generates a matrix of size p x p with ones on
  ##the diagInd diagonal. Restriction is abs(diagInd) <= p. Returns 0 matrix otw.
  
  diagNew <- abs(diagInd);
  if(diagNew <= p){
    firstPart <- diag(1,p,p-diagNew);
    secondPart <- array(0,dim=c(p,diagNew));
    
    A <- cbind(secondPart,firstPart);
    if(diagInd < 0){
      A <- t(A);
    }
  }
  else{
    A <- c(0)
  }
  A
}
## end diagExtend



