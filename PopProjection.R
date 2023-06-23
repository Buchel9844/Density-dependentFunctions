# Matlab:
# dndt = R.*n.*(1-exp(- sum(a.*n_array,2)  )) - (b*n).*n;

library(deSolve)

Ricker_function1 <- function(time,
                  state,
                  pars) {
  
  S <- length(pars$r) # number of species
  n <- state[1:S] # species densities
  r <- pars$r # maximum growth rates 
  a <- pars$a # beneficial interactions
  b <- pars$b # harmful interactions 
  
  aii <- diag(a)
  aij <- a
  diag(aij) <- 0
  
  aii <- 0.2
  aij <- Amin
  
  Fi <- lambda* exp(aii%*%Ni + aij%*%Nj)
  Fj <- lambda* exp(ajj%*%Nj + aji%*%Ni)
  
  dNidt <- (1-g)*s + g* Fi
  dNjdt <- (1-g)*s + g* Fj
  
  dndt <- r*n*( 1 - exp(- aii - aij%*%n) ) - n*(b%*%n) # equation for change in abundances
  
  return(list(dndt))
}

pars <- list(r = c(2.3, 2.3),
             a = matrix(c(0.05, 0.1, 0.05, 0.3), nrow = 2, byrow = T),
             b = matrix(c(0.1, 0.05, 0.2, 0.05), nrow = 2, byrow = T))

state <- c(10, 10) # initial species densities

time = seq(0, 100, by = 1)

ode(y = state, times = time, func = baker, parms = pars)

