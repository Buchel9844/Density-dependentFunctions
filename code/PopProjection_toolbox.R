alpha_function2 <- function(Amin, Aslopes,N,N0){
  alpha = Amin + Aslopes*(N-N0)
  if((N-N0) >10 ){
    alpha = Amin + Aslopes*(10)
  }
  return(alpha)
}
alpha_function3 <- function(Amin, Aslopes,c,N,N0){
  alpha = Amin + Aslopes*(1-exp(-c*(N-N0)))
  if((N-N0) >10 ){
    alpha = Amin + Aslopes*(1-exp(-c*(10)))
  }
  return(alpha)
}
alpha_function4  <- function(Amin, Aslopes,c,N,N0){
  e = exp(-Aslopes*(N-N0)) # c is stretching the graph horizontally 
  if((N-N0) >10 ){
    e = exp(-Aslopes*(10))
  }
  a = c*(1-e) #stretching the graph vertically
  d = Amin
  alpha = (a/(1 + e)) + d
  
  return(alpha)
}
library(deSolve)
alpha_function4(-0.6380506,0.1211577,0.8000038,5,1)

Ricker_function <- function(t,
                            y,
                            p) {
  
  Ni <- y[1] # species i initial densities
  Nj <- y[2] # species j initial densities
  with( as.list(p),{
  
  if(function.int==1){
    aii <- a_initial[1,1]
    aij <- a_initial[1,2]
    aji <- a_initial[2,1]
    ajj <- a_initial[2,2]
  }
  if(function.int==2){
    aii <- alpha_function2(a_initial[1,1], a_slope[1,1],g[1]*Ni, Nmax[1,1])
    aij <- alpha_function2(a_initial[1,2], a_slope[1,2],g[2]*Nj, Nmax[1,2])
    aji <- alpha_function2(a_initial[2,1], a_slope[2,1],g[1]*Ni, Nmax[2,1])
    ajj <- alpha_function2(a_initial[2,2], a_slope[2,2],g[2]*Nj, Nmax[2,2])
  }
  if(function.int==3){
    aii <- alpha_function3(a_initial[1,1], a_slope[1,1],c[1,1],g[1]*Ni, Nmax[1])
    aij <- alpha_function3(a_initial[1,2], a_slope[1,2],c[1,2],g[2]*Nj, Nmax[2])
    aji <- alpha_function3(a_initial[2,1], a_slope[2,1],c[2,1],g[1]*Ni, Nmax[1])
    ajj <- alpha_function3(a_initial[2,2], a_slope[2,2],c[2,2],g[2]*Nj, Nmax[2])
  }
  if(function.int==4){
    aii <- alpha_function4(a_initial[1,1], a_slope[1,1],c[1,1],g[1]*Ni, Nmax[1])
    aij <- alpha_function4(a_initial[1,2], a_slope[1,2],c[1,2],g[2]*Nj, Nmax[2])
    aji <- alpha_function4(a_initial[2,1], a_slope[2,1],c[2,1],g[1]*Ni, Nmax[1])
    ajj <- alpha_function4(a_initial[2,2], a_slope[2,2],c[2,2],g[2]*Nj, Nmax[2])
  }
  
  
  Fi <-  exp(lambda[1] + aii * g[1] * Ni + aij * g[2] * Nj)
  Fj <-  exp(lambda[2] + ajj * g[2] * Nj + aji * g[1] * Ni)
  
  dNi <- (1-g[1]) * s[1] + g[1] * Fi
  dNj <- (1-g[2]) * s[2] + g[2] * Fj
  
  return(list(c(dNi,dNj))) })
}

Ricker_solution <- function(gens,
                            state,
                            pars) {
  Nmax <- pars$Nmax # density at which fecundity is max - effect of neighbors is 0
  g <- pars$g # germination rate 
  s <- pars$s #seed survival
  lambda <- pars$lambda # intrinsic growth rate
  function.int <- pars$function.int # which int.function
  a_initial <- pars$a_initial # which int.function
  a_slope <- pars$a_slope # which int.function
  c <- pars$c # which int.function

  df <- data.frame( t=0:gens,  Ni=numeric(1+gens),  Nj =numeric(1+gens) )
  df[1,2:3] <- c(state[1],state[2]) #species i initial densities
  
  for(t in 1:gens){

  Ni <- df[t,"Ni"] # species i densities
  Nj <- df[t,"Nj"] # species j  densities
  
  if(function.int==1){
    aii <- a_initial[1,1]
    aij <- a_initial[1,2]
    aji <- a_initial[2,1]
    ajj <- a_initial[2,2]
  }
  if(function.int==2){
    aii <- alpha_function2(a_initial[1,1], a_slope[1,1],g[1]*Ni, Nmax[1,1])
    aij <- alpha_function2(a_initial[1,2], a_slope[1,2],g[2]*Nj, Nmax[1,2])
    aji <- alpha_function2(a_initial[2,1], a_slope[2,1],g[1]*Ni, Nmax[2,1])
    ajj <- alpha_function2(a_initial[2,2], a_slope[2,2],g[2]*Nj, Nmax[2,2])
  }
  if(function.int==3){
    aii <- alpha_function3(a_initial[1,1], a_slope[1,1],c[1,1],g[1]*Ni, Nmax[1])
    aij <- alpha_function3(a_initial[1,2], a_slope[1,2],c[1,2],g[2]*Nj, Nmax[2])
    aji <- alpha_function3(a_initial[2,1], a_slope[2,1],c[2,1],g[1]*Ni, Nmax[1])
    ajj <- alpha_function3(a_initial[2,2], a_slope[2,2],c[2,2],g[2]*Nj, Nmax[2])
  }
  if(function.int==4){
    aii <- alpha_function4(a_initial[1,1], a_slope[1,1],c[1,1],g[1]*Ni, Nmax[1])
    aij <- alpha_function4(a_initial[1,2], a_slope[1,2],c[1,2],g[2]*Nj, Nmax[2])
    aji <- alpha_function4(a_initial[2,1], a_slope[2,1],c[2,1],g[1]*Ni, Nmax[1])
    ajj <- alpha_function4(a_initial[2,2], a_slope[2,2],c[2,2],g[2]*Nj, Nmax[2])
  }

  
  Fi <-  exp(lambda[1] + aii * g[1]*Ni + aij *g[2]*Nj)
  Fj <-  exp(lambda[2] + ajj * g[2]*Nj + aji *g[1]*Ni)
  
  Nit1 <- ((1-g[1]) * s[1] + g[1] * Fi)*Ni
  Njt1<- ((1-g[2]) * s[2] + g[2] * Fj)*Nj
  df[t+1,2:3] <- c(Nit1, Njt1)
  }
  return(df)
}
