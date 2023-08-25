library(forecast)
library(tseries)
library(plotly)
##########################################################################################################
# Stability metric
##########################################################################################################
 df.stability <- df.sim[which(df.sim$invader =="both"),]
 
 summary.df.stability <- NULL
 nsims <- 1000
 t.num = 100 # number of generation
 
 for(i in 1:nsims){
   for( function.int in 1:4){
     print(paste0("int ", i,"for funct ",function.int))
   df.stability.n <-  df.stability[which(df.sim$sim == i &
                                       df.sim$function.int == function.int)]
 
   df.stability.n <- df.stability.n[c(10:t.num*100),] # burn first 10 generations
 
 
 #Mean Variance over time
 # Plaza et al,2012 https://doi.org/10.1111/j.1654-1103.2011.01381.x
 #  amplitude of population fluctuations by means of the standard deviation
      msdi <- sd(log10(df.stability.n$Ni))
      msdj <- sd(log10(df.stability.n$Nj))
      
#ARIMA is the abbreviation for AutoRegressive Integrated Moving Average.      
      fitARIMA <- auto.arima(log(test$Nj)) # function in R uses a combination of unit root tests, minimization of the AIC and MLE to obtain an ARIMA model
      coeftest(fitARIMA)
      
    df.stability.summary.n <- data.frame(sim == i,
                                     function.int == function.int,
                                     species = c("Ni","Nj"),
                                     mean=c(mean(df.sum.small$Ni),mean(df.sum.small$Nj)),
                                     msd = c(msdi,msdj))
    
    df.stability.summary <- bind_rows(df.stability.summary,df.stability.summary.n)
                                     
 
   }
}

stand.variable <- c("mean","median","var","upperbound","lowerbound")


summary.df.abundance  %>%
  gather(mean,median,var,upperbound,lowerbound ,key="summary",value="value") %>%
  ggplot(aes(x=as.factor(species), y= value)) + 
  geom_boxplot(aes(fill = as.factor(function.int)))+
  facet_grid(summary~., scales = "free") + theme_bw() + 
  scale_fill_colorblind() + ylim(0,10)


ggsave(last_plot(),
       file = "figures/summary.abundances.pdf")


#---- Visualise for one sim and one int.function the different graph associated with the above metrics
i = 2
function.int = 4

test <- df.stability[which(df.sim$sim == i &
                             df.sim$function.int == function.int)]
testsmall <- as.matrix(test[c(10:t.num),c("Ni","Nj")])


# Turchin and Taylor, 1992
# compute AFC : the correlation between numbers at one point in time and those at different times in the past
# negative curve: a stationary process with exponential return to equilibrium
# negative exp curve in the negative :  process with nonstationary mean and no periodicity
# constant waves : a process driven by an exogenous periodic force, or phase-remembering quasi-cycle
# decreasing waves: a stationary process with endogenously generated periodicity, or phase-forgetting quasi-cycle
acf(log10(test$Nj)) 
pacf(log10(test$Nj))
acf(log10(test$Ni))
pacf(log10(test$Ni))

#ARIMA is the abbreviation for AutoRegressive Integrated Moving Average. 
# Auto Regressive (AR) terms refer to the lags of the differenced series: p = the number of autoregressive term 
# Moving Average (MA) terms refer to the lags of errors: q 
# I is the number of difference used to make the time series stationary.
#d = the number of nonseasonal differences
#q = the number of moving-average terms
#ARIMA(p,d,q)”
# interpret results
#  In a model with 2 or more AR coefficients, the sum of the coefficients determines the speed of mean reversion, and the series may also show an oscillatory pattern
adf.test(test$Ni)

fitARIMA <- auto.arima(test$Nj) # function in R uses a combination of unit root tests, minimization of the AIC and MLE to obtain an ARIMA model
coeftest(fitARIMA)
summary(fitARIMA)

# visualise that the remaining values are random noise
acf(fitARIMA$residuals)
qqnorm(fitARIMA$residuals)
qqline(fitARIMA$residuals)

pcor(testsmall,method = c("pearson")) # high value shows that hat x and y are highly consistent and they increase with each other

# Plaza et al,2012 https://doi.org/10.1111/j.1654-1103.2011.01381.x
#  amplitude of population fluctuations by means of the standard deviation
sd(log10(test$Nj))
sd(log10(test$Ni))

#Berryman and Lima, 2007 https://doi.org/10.1890/06-0609.1. - limitation of AFC - check following graph conitnuity
plot(test$dNi,c(test$Ni[-1],NA))
plot(test$dNj,c(test$Nj[-1],NA))


library(plotly)
p <- plot_ly(z = volcano, type = "surface")
p 
vreq_classic(data.frame(x=acf(test$Nj),y=acf(test$Nj)))

# Zhao and al,2020 Package https://cran.r-project.org/web/packages/tsvr/vignettes/tsvrvignette.pdf
# > 1 synchrony | < 1 asynchrony
testsmall <- t(as.matrix(test[c(10:t.num),c("Ni","Nj")]))


vreq_classic(testsmall)
vr.trial <- tsvreq_classic(testsmall)
summary(vr.trial)
plot(vr.trial)
aggresShort <- aggts(vr.trial, vr.trial$ts[vr.trial$ts<4])
aggresLong <- aggts(vr.trial, vr.trial$ts[vr.trial$ts>=4])
vr.trial.short <- aggresShort[[3]]
vr.trial.long <- aggresLong[[3]]
temp <-  vreq_classic(testsmall)
vr.class<- temp[[3]]
plot(c(1,2,3), c(vr.trial.short, vr.trial.long, vr.class), 
     col=c("blue", "darkgreen", "lightseagreen"), pch=16, cex=2, xaxt="n", yaxt="n", ylim=c(-.05, 2.05), xlim=c(0.7,3.3))
axis(side=1, at=c(1,2,3), labels=c("Short \n Timescale", "Long \n Timescale", "Classic \n"), tick=FALSE, cex=1.25)
axis(side=2, at=c(0,1,2), labels=TRUE, las=1, tick=TRUE, hadj=-.25, lwd.ticks = 1, tck=-.03)
mtext("Variance Ratio", side=3, outer=FALSE, line=1)
abline(h=1, col="grey", lty=2, lwd=2)
text(x=0.72, y= 1.98, "g)", cex=1.25)


cov(testsmall$Nj,testsmall$Ni)/var(testsmall$Ni) 
cov(testsmall$Ni,testsmall$Nj)/var(testsmall$Nj)

cov(testsmall$total,testsmall$Ni)/var(testsmall$Ni)
cov(testsmall$total,testsmall$Nj)/var(testsmall$Nj)

##########################################################################################################
# CHAOS metric
##########################################################################################################
install.packages(mgcv)
library(mgcv) # for GAM ? 
#https://cran.r-project.org/web/packages/mgcv/mgcv.pdf


library(tseriesCHAOS) # for lyapunov state 
#https://cran.r-project.org/web/packages/nonlinearTseries/nonlinearTseries.pdf
