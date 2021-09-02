library(spatstat)
library(R.matlab)
library(RandomFields)
library(pracma)

# generating different environments by varying bacterial dispersion (var_p)
# average intensity (number of cells) = exp(mu + var/2) (as var changes, mu needs to be adjusted 
# to keep intensity fixed)
var_p <- 1.9
mu_p <- -7.52
scale_p <- 25
xwin = c(0,1000)
ywin = c(0,3000)
vec = c(logspace(-3,1.5,10),80)
target_val = exp(mu_p + var_p/2)*3001*1001 # expected number of bacteria in sim domain
nsim = 2000

for(val in vec)
{
  var_p <- 1.9 - (1.9-val)
  mu_p <- -7.52 + (1.9-val)/2
  X <- rLGCP("exp", mu=mu_p, var=var_p, scale=scale_p, win=owin(xwin, ywin), nsim=nsim)
  total = c()
  for (ii in 1:nsim)
  {
    ind <- paste("Simulation",ii)
    total <- c(total,X[[ind]][["n"]])
  }
  # pick the simulation whose bacteria count is closest to expectation
  X <- X[[paste("Simulation",which.min(abs(total-target_val)))]]
  print(X[["n"]])
  plot(X)
  soil_LGCP <- as.matrix(as.data.frame(X))
  writeMat(paste("soil_var_",round(val,3),".mat",sep=""), soil_LGCP=soil_LGCP, var_p=var_p, 
           mu_p=mu_p, scale_p=scale_p, xwin=xwin, ywin=ywin)
}