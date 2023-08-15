alpha_function2 <- function(Amin, Aslopes,N,N0){
  alpha = Amin + Aslopes*(N-N0)
  #if((N-N0) >10 ){
  #  alpha = Amin + Aslopes*(10)
  #}
  return(alpha)
}
alpha_function3 <- function(Amin, Aslopes,c,N,N0){
  alpha = Amin + Aslopes*(1-exp(-c*(N-N0)))
  #if((N-N0) >10 ){
  #  alpha = Amin + Aslopes*(1-exp(-c*(10)))
  #}
  return(alpha)
}
alpha_function4  <- function(Amin, Aslopes,c,N,N0){
  e = exp(-c*(N-N0)) # c is stretching the graph horizontally 
  #if((N-N0) >10 ){
  #  e = exp(-Aslopes*(10))
  #}
  a = Aslopes*(1-e) #stretching the graph vertically
  d = Amin
  alpha = (a/(1 + e)) + d
  
  return(alpha)
}
Ricker_solution_ODE <- function(gens,
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
  
  df <- data.frame( t=0:gens,  Ni=numeric(1+gens),  Nj =numeric(1+gens) ,
                    dNi=numeric(1+gens),  dNj =numeric(1+gens) )
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
    Njt1 <- ((1-g[2]) * s[2] + g[2] * Fj)*Nj
    Nidt <- Nit1/Ni
    Njdt <- Njt1/Nj
    
    df[t+1,2:5] <- c(Nit1, Njt1, Nidt, Njdt)
  }
  names(df) <- c("time","Ni","Nj","dNi","dNj")
  return(df)
}


Ricker_solution <- function(gens,
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
    Njt1 <- ((1-g[2]) * s[2] + g[2] * Fj)*Nj
    df[t+1,2:3] <- c(Nit1, Njt1)
  }
  names(df) <- c("time","Ni","Nj")
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


# function modified from https://github.com/laurenmh/avena-erodium/blob/master/invader_resident_comparison.R
grwr = function(par.dat, t.num,function.int) {
  
  N1 = c(1,0) #c("Ni", "Nj")
  Niequil <- Ricker_solution_mono(state =  N1, pars= par.dat, gens=t.num,
                                  function.int)
  
  N2 = c(0,1) #c("Ni", "Nj")
  Njequil <- Ricker_solution_mono(state =  N2, pars= par.dat, gens=t.num,
                                  function.int)
  
  
  N3 = c(Niequil[t.num, 2],1) #c("Ni*", "Nj")
  Njinvade <- Ricker_solution(state = N3, pars= par.dat, gens=t.num,
                              function.int)
  
  Njinvade$invader <- "Nj"
  Njinvade$time <- as.numeric(row.names(Njinvade))
  
  N4 = c(1, Njequil[t.num, 3]) #c("Ni", "Nj")
  Niinvade <- Ricker_solution(state = N4, pars= par.dat, gens=t.num,
                              function.int)
  
  Niinvade$invader <- "Ni"
  Niinvade$time <- as.numeric(row.names(Niinvade))
  
  Niinvade$grwr <- NA
  Niinvade$Cgrwc<- NA
  Njinvade$grwr<- NA
  Njinvade$Cgrwc<- NA
  Niinvade$grwrChesson<- NA
  Njinvade$grwrChesson<- NA
  Niinvade$grwrStouffer<- NA
  Njinvade$grwrStouffer<- NA
  
  for( t in 1:t.num){
    Niinvade$grwr[t] <- Niinvade$Ni[t+1]/Niinvade$Ni[t]
    Niinvade$Cgrwc[t] <- Niinvade$Nj[t+1]/Niinvade$Nj[t]
    Njinvade$grwr[t] <- Njinvade$Nj[t+1]/Njinvade$Nj[t]
    Njinvade$Cgrwc[t] <- Njinvade$Ni[t+1]/Njinvade$Ni[t]
    
    Niinvade$grwrChesson[t] = log(Niinvade$grwr[t]) - log(Njinvade$Cgrwc[t])# for Ni
    Njinvade$grwrChesson[t] = log(Njinvade$grwr[t]) - log(Niinvade$Cgrwc[t]) # for Nj
    
    Niinvade$grwrStouffer[t] = log(Niinvade$grwr[t]) 
    Njinvade$grwrStouffer[t] = log(Njinvade$grwr[t]) 
  }
  
  vert <- rbind(Niinvade, Njinvade) 

  return(list(vert=vert))
}  
