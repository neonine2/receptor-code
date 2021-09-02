library(spatstat)
library(R.matlab)
library(RandomFields)
library(pracma)

# homogeneous LGCP with exponential covariance function
# values reported in Raynaud and Nunan for 10^9 cells per g of soil
var_p <- 1.9
mu_p <- -7.52
scale_p <- 25
xwin <- c(0,1000)
ywin <- c(0,3000)
target_val <- exp(mu_p + var_p/2)*3001*1001 # expected number of bacteria in sim domain
nsim <- 2000

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
# writing all bacteria position into Matlab file
writeMat("soil_var_2.mat", soil_LGCP=soil_LGCP, var_p=var_p, 
         mu_p=mu_p, scale_p=scale_p, xwin=xwin, ywin=ywin)
