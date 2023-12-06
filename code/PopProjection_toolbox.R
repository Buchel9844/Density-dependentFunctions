alpha_function2 <- function(Amin, Aslopes,N,N0){
  alpha = Amin + Aslopes*(N-N0)
  #if((N-N0) >10 ){
  #  alpha = Amin + Aslopes*(10)
  #}
  return(alpha)
}
alpha_function3 <- function(Amin, Aslopes,c,N,N0){
  alpha = Amin + c*(1-exp(Aslopes*(N-N0)))
  #if((N-N0) >10 ){
  #  alpha = Amin + Aslopes*(1-exp(-c*(10)))
  #}
  return(alpha)
}
alpha_function4  <- function(Amin, Aslopes,c,N,N0){
  e = exp(Aslopes*(N-N0)) # c is stretching the graph horizontally 
  #if((N-N0) >10 ){
  #  e = exp(-Aslopes*(10))
  #}
  a = c*(1-e) #stretching the graph vertically
  d = Amin
  alpha = (a/(1 + e)) + d
  
  return(alpha)
}
Ricker_solution_withalpha <- function(gens,
                                     state,
                                     pars,
                                     function.int,
                                     add_external_factor) {
  Nmax <- pars$Nmax # density at which fecundity is max - effect of neighbors is 0
  g <- pars$g # germination rate 
  s <- pars$s #seed survival
  lambda <- pars$lambda # intrinsic growth rate
  a_initial <- pars$a_initial # which int.function
  a_slope <- pars$a_slope # which int.function
  c <- pars$c # which int.function
  
  df <- data.frame( t=0:gens,  Ni=numeric(1+gens),  Nj =numeric(1+gens) ,
                    dNi=numeric(1+gens),  dNj =numeric(1+gens) )
  df[1,2:3] <- c(state[1],state[2]) #species i initial densities
  
  if(add_external_factor =="none"){
    e <- rep(0,times=gens)
  }
  if(add_external_factor =="season"){
    e <- pars$e_seasonal
  }
  if(add_external_factor =="noise"){
    e <- pars$e_noise
  }
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
    
    
    Fi <-  exp(lambda[1] + aii * g[1]*Ni + aij *g[2]*Nj + e[t])
    Fj <-  exp(lambda[2] + ajj * g[2]*Nj + aji *g[1]*Ni + e[t])
    
    Nit1 <- ((1-g[1]) * s[1] + g[1] * Fi)*Ni 
    Njt1 <- ((1-g[2]) * s[2] + g[2] * Fj)*Nj
    Nidt <- Nit1/Ni
    Njdt <- Njt1/Nj
    
    df[t,4:5] <- c(Nidt, Njdt)
    df[t,6:9] <- c(aii,ajj,aij,aji)
    df[t+1,2:3] <- c(Nit1, Njt1)
  }
  df[c(1:length(e)),10] <- c(e)
  names(df) <- c("time","Ni","Nj","dNi","dNj","aii","ajj","aij","aji","external.fact")
  return(df)
}


Ricker_solution <- function(gens,
                                state,
                                pars,
                           function.int,
                           add_external_factor) {
  Nmax <- pars$Nmax # density at which fecundity is max - effect of neighbors is 0
  g <- pars$g # germination rate 
  s <- pars$s #seed survival
  lambda <- pars$lambda # intrinsic growth rate
  a_initial <- pars$a_initial # which int.function
  a_slope <- pars$a_slope # which int.function
  c <- pars$c # which int.function
  
  df <- data.frame( t=0:gens,  Ni=numeric(1+gens),  Nj =numeric(1+gens) ,
                    dNi=numeric(1+gens),  dNj =numeric(1+gens) )
  df[1,2:3] <- c(state[1],state[2]) #species i initial densities
  
  if(add_external_factor =="none"){
    e <- rep(0,times=gens)
  }
  if(add_external_factor =="season"){
    e <- pars$e_seasonal
  }
  if(add_external_factor =="noise"){
    e <- pars$e_noise
  }
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

    
    Fi <-  exp(lambda[1] + aii * g[1]*Ni + aij *g[2]*Nj + e[t])
    Fj <-  exp(lambda[2] + ajj * g[2]*Nj + aji *g[1]*Ni + e[t])
    
    Nit1 <- ((1-g[1]) * s[1] + g[1] * Fi)*Ni 
    Njt1 <- ((1-g[2]) * s[2] + g[2] * Fj)*Nj
    Nidt <- Nit1/Ni
    Njdt <- Njt1/Nj
    
    df[t,4:5] <- c(Nidt, Njdt)
    df[t+1,2:3] <- c(Nit1, Njt1)
  }
  df[c(1:length(e)),6] <- c(e)
  names(df) <- c("time","Ni","Nj","dNi","dNj","external.fact")
  return(df)
}

Ricker_solution_mono <- function(gens,
                                 state,
                                 pars,
                                 function.int) {
  Nmax <- pars$Nmax # density at which fecundity is max - effect of neighbors is 0
  g <- pars$g # germination rate 
  s <- pars$s #seed survival
  lambda <- pars$lambda # intrinsic growth rate
  a_initial <- pars$a_initial # which int.function
  a_slope <- pars$a_slope # which int.function
  c <- pars$c # which int.function
  
  if(state[2]==0){ # monoculture of Ni
    df <- data.frame( t=0:gens,  Ni=numeric(1+gens))
    
    df[1,2] <- c(state[1]) #species i initial densities
    
    for(t in 1:gens){
      
      Ni <- df[t,"Ni"] # species i densities
      
      if(function.int==1){
        aii <- a_initial[1,1]
        
      }
      if(function.int==2){
        aii <- alpha_function2(a_initial[1,1], a_slope[1,1],g[1]*Ni, Nmax[1,1])
      }
      if(function.int==3){
        aii <- alpha_function3(a_initial[1,1], a_slope[1,1],c[1,1],g[1]*Ni, Nmax[1])
      }
      if(function.int==4){
        aii <- alpha_function4(a_initial[1,1], a_slope[1,1],c[1,1],g[1]*Ni, Nmax[1])
      }
      
      Fi <-  exp(lambda[1] + aii * g[1]*Ni )
      
      Nit1 <- ((1-g[1]) * s[1] + g[1] * Fi)*Ni
      
      df[t+1,2] <- c(Nit1)
    }
    names(df) <- c("time","Ni")
  }
  if(state[1]==0){ # monoculture of Nj
    df <- data.frame( t=0:gens,  Nj=numeric(1+gens))
    
    df[1,2] <- c(state[2]) #species i initial densities
    
    for(t in 1:gens){
      Nj <- df[t,"Nj"] # species j  densities
      
      if(function.int==1){
        
        ajj <- a_initial[2,2]
      }
      if(function.int==2){
        
        ajj <- alpha_function2(a_initial[2,2], a_slope[2,2],g[2]*Nj, Nmax[2,2])
      }
      if(function.int==3){
        
        ajj <- alpha_function3(a_initial[2,2], a_slope[2,2],c[2,2],g[2]*Nj, Nmax[2])
      }
      if(function.int==4){
        ajj <- alpha_function4(a_initial[2,2], a_slope[2,2],c[2,2],g[2]*Nj, Nmax[2])
      }
      
      Fj <-  exp(lambda[2] + ajj * g[2]*Nj)
      
      Njt1<- ((1-g[2]) * s[2] + g[2] * Fj)*Nj
      
      df[t+1,2:3] <- c( Njt1)
    }
    names(df) <- c("time","Nj")
  }
  
  return(df)
}


Ricker_solution_NatData <- function(gens,
                                    state,
                                    pars) {
  function.int <- pars$function.int[1]
  Nmax <- pars$Nmax # density at which fecundity is max - effect of neighbors is 0
  g <- pars$g[1] # germination rate 
  s <- pars$s[1] #seed survival
  lambda <- pars$lambda[1] # intrinsic growth rate
  a_initial <- pars$alpha_init # which int.function
  a_slope <- pars$alpha_slope # which int.function
  c <- pars$alpha_c # which int.function
   
  df <- state
  
  for(t in 1:(gens-1)){
    Nt1 <- c()
    Fec <- c()
    Nt <-  df[t,]  # species i densities
    if(function.int==1){
        a <- a_initial
      }
      if(function.int==2){
        a <- alpha_function2(a_initial, a_slope,g*Nt, Nmax)
      }
      if(function.int==3){
        a <- alpha_function3(a_initial, a_slope,c,g*Nt, Nmax)
      }
      if(function.int==4){
        a <- alpha_function4(a_initial, a_slope,c,g*Nt, Nmax)
      }
      
      Fec <- exp(lambda + sum(a*g*Nt))
      
      Nt1 <- ((1-g) * s + g* Fec)*Nt$conspecific
      df$conspecific[t+1] <- c(Nt1)
    }
    
  df$time <- c(1:gens) 
  return(df)
}
# function modified from https://github.com/laurenmh/avena-erodium/blob/master/invader_resident_comparison.R

GrowthSimInv = function(par.dat, t.num,function.int,add_external_factor) {
  
  N1 = c(1,0) #c("Ni", "Nj")
  Niequil <- Ricker_solution_mono(state =  N1, pars= par.dat, gens=t.num,
                                  function.int)
  
  N2 = c(0,1) #c("Ni", "Nj")
  Njequil <- Ricker_solution_mono(state =  N2, pars= par.dat, gens=t.num,
                                  function.int)
  
  
  N3 = c(Niequil[t.num, 2],1) #c("Ni*", "Nj")
  Njinvade <- Ricker_solution(state = N3, pars= par.dat, gens=t.num,
                              function.int,add_external_factor)
  
  Njinvade$invader <- "Nj"
  Njinvade$time <- as.numeric(row.names(Njinvade))
  
  N4 = c(1, Njequil[t.num, 3]) #c("Ni", "Nj")
  Niinvade <- Ricker_solution(state = N4, pars= par.dat, gens=t.num,
                              function.int,add_external_factor)
  
  Niinvade$invader <- "Ni"
  Niinvade$time <- as.numeric(row.names(Niinvade))
  
  df <- rbind(Niinvade, Njinvade) 
  
  return(df)
}  



