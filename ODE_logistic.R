# from https://hankstevens.github.io/Primer-of-Ecology/DDgrowth.html


clogistic <- function(time, y, parameters){
  # y is a single value of N; if we have more than one population, 
  # e.g. competition, then y would be a vector with two elements, one for 
  # each population.
  
  # parameters is a vector of...parameters....
  
  # Both y and parameters have names for each element,
  # e.g., y <- c(N=10) - see the text.
  
  with(as.list(c(y, parameters)), {
    # "with" creates an "environment" within which R looks for values
    # c() combines the two vectors into one
    # as.list() turns the resulting vector into a list
    
    # our differential equation for growth
    dN.dt <- r * N * (1 - alpha * N)
    
    # making the function create output, which must be a list with
    # one element
    return( list( dN.dt ) )
  } )
}


logistic2 <- function(t, y, p){
  dN.dt <- p[1] * y[1] * (1 - p[2] * y[1])
  return( list( dN.dt ) }
          
p <- c(r=1, alpha = 0.2)
y0 <- c(N=2)
t <- 0:20        

out <- ode(y=y0, times=t, func=clogistic, parms=p)
out[1:5,]
plot(out)
