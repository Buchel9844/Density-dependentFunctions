##########################################################################################################
# Stability metric
##########################################################################################################


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
#  In a model with 2 or more AR coefficients, the sum of the coefficients determines the speed of mean reversion, 
# and the series may also show an oscillatory pattern
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

##########################################################################################################
# CHAOS metric. for continuous time
##########################################################################################################
#----- Define jacobian with ODE - UNSUCCESSFUL ---- 

p[["function.int"]] <- 4 # function.int
state= c(Ni=1,Nj=1)
N.init <- c(1,5,10)
out <- NULL
for(i in  1:3){
  # project the population
  state <- c(Ni=N.init[i],Nj=N.init[i])
  
  # project the population
  a <- deSolve::lsode(y = state, times = seq(1,t.num*10,by=1), 
                      func = Ricker_solution_ODE, parms = p, method = "euler")
  
  b <- data.frame(a, N0=as.factor(N.init[i])) # b now has cols: t, N, rd
  # store data in rows of out "bind rows" with rbind
  out <- rbind.data.frame(out, b)
}
test_plot <-  out %>% 
  gather(Ni,Nj, key="pop", value="abundance") %>% 
  ggplot(aes(x = time, y = abundance)) +
  geom_line(aes(colour = pop, linetype=N0)) +
  #geom_point(aes(color = pop, shape=N0)) +
  scale_color_brewer(NULL, palette = "Set1") +
  labs(title = "Ricker model",
       #subtitle = paste(names(pars), pars, sep = " = ", collapse = "; "),
       x = "Time", y = "Population density")
test_plot

n <- out[which(out$N0 == 1), 2]
n <- n[901:1000]
ntp1 <- n[-1]
nt <- n[-length(n)]
qplot(nt, ntp1, geom=c("point", "path") ) + 
  labs(y="N[t+1]", x="N[t]")


#----- Define Lyapunov exponent with DATA ---- 
# with FORECAST 
library(forecast)
findfrequency(log(df.stability.n$Ni))
decomp(log(df.stability.n$Ni))

x <- log(df.stability.n$Ni)

N <- length(x)
freq <- findfrequency(x)
fx <- c(frequency=(exp((freq-1)/50)-1)/(1+exp((freq-1)/50)))
x <- ts(x,f=freq)

if(freq > N-10)
  stop("Insufficient data")
Ly <- numeric(N-freq)
for(i in 1:(N-freq)){ # function from https://robjhyndman.com/hyndsight/tscharacteristics/
  idx <- order(abs(x[i] - x))
  idx <- idx[idx < (N-freq)]
  j <- idx[2]
  Ly[i] <- log(abs((x[i+freq] - x[j+freq])/(x[i]-x[j])))/freq
  if(is.na(Ly[i]) | Ly[i]==Inf | Ly[i]==-Inf)
    Ly[i] <- NA
}
Lyap <- mean(Ly,na.rm=TRUE)
fLyap <- exp(Lyap)/(1+exp(Lyap))


# with GAM Generalised additive models - 
install.packages(mgcv)
library(mgcv) # for GAM ? 
#https://cran.r-project.org/web/packages/mgcv/mgcv.pdf
df.stability.n.gam <- df.stability.n %>%
  mutate(y= log(Ni),
         x = c(log(Ni[-1]),NA),
         Njt2 = c(Nj[-1],NA))

gam.test <- gam( y ~ s(x), data = df.stability.n.gam,
                 method = "REML")
plot(gam.test)

model_matrix <- predict(gam.test, type = "lpmatrix")
plot(y~x,data = df.stability.n.gam)
abline(h = 0)
lines(x[-1], model_matrix[, "s(x).1"], type = "l", lty = 2)
lines(x[-1], model_matrix[, "s(x).2"], type = "l", lty = 2)


x_new <- seq(0, 100, length.out = 100)
y_pred <- predict(gam.test, data.frame(time = x_new))

ggplot(df.stability.n[c(1:100),], aes(y=Ni, x=time)) +
  geom_point() +
  geom_smooth(method = "gam",formula = y ~ s(x))
par(mfrow = c(2, 2))
gam.check(gam.test)

# with DChaos: Chaotic Time Series Analysis- 

install.packages("DChaos")
library(DChaos)

df.stability.n <-  df.stability[which(df.stability$sim== 2 &
                                        df.stability$function.int == 4
                                      & df.stability$time > 10),]

df.stability.n %>%
  gather(Ni,Nj, key="pop", value="abundance") %>% 
  ggplot(aes(x = time, y = log(abundance))) +
  geom_line(aes(colour = pop)) +
  #geom_point(aes(color = pop, shape=N0)) +
  scale_color_brewer(NULL, palette = "Set1") +
  theme_bw()

Exponent <- DChaos::lyapunov(log(df.stability.n$Ni),timelapse="FIXED")
summary(Exponent)
Lyap <- mean(Exponent$exponent.mean[,"Estimate"])
Lyap
fLyap <- exp(Lyap)/(1+exp(Lyap))
fLyap