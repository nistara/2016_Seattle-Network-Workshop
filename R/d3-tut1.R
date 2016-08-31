
##
## Network Modeling for Epidemics 2016
## Day 3 | Tutorial 1: SIS Epidemic in a One-Mode Network
##

## Load EpiModel
library(EpiModel)


# 1. Model Estimation -----------------------------------------------------
## -----------------------------------------------------------------------------

## Initialize the network
nw <- network.initialize(n = 500, directed = FALSE) # undirected coz we're simulating
# sexual contact

## Formation model formula
formation <- ~edges + concurrent
# term for edges (intercept - density of network)
# how many people have degree 2 or greater - concurrent


## Target statistics for formation model
target.stats <- c(175, 110)
# correspondng to mean degree of 0.7
# 22% edges exhibit concurrency (i.e. have 2 or greater degrees)

## Constraints on the degree distribution
constraints <- ~bd(maxout = 3)
# bd = by degree, where the maximum is 3

## Calculate the dissolution coefficients
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)
coef.diss # average duration of partnerships.
# log od duration - 1 = 3.8918
# use offset because they're fixed?


## Show arguments to netest function
args(netest)

## Estimate the temporal ERGM
est <- netest(nw, formation, target.stats, coef.diss, constraints)

## Run diagnostics on the model fit
dx <- netdx(est, nsims = 10, nsteps = 1000,
            nwstats.formula = ~edges + meandeg + degree(0:4) + concurrent)
# Note, we have additional terms in this mode, to see how they're doing
# meandeg - mean degree
# degree(0:4) - no of nodes with 0, 1, 2, 3, 4 ties
# Putting degree 4 is actually a check - coz we had a constraint
# Formation diag - summaru statistics of the cross section of network
# 10,000 observations of networks in total, and the summary is the mean of these
# 10*1000 = 10000


## Print the diagnostics
dx

## Plot the diagnostics for the formation model
plot(dx)
# the dashed line - target statistic
# Band - interquartile range
?plot.netdx

## Plot the dissolution model diagnostics
par(mfrow = c(1, 2))
plot(dx, type = "duration", mean.col = "black")
plot(dx, type = "dissolution", qnts = 0.5,
     mean.lines = FALSE, sim.lines = FALSE)

## An example of static diagnostics
dx.static <- netdx(est, nsims = 10000, dynamic = FALSE)
dx.static
# if average duration of edges is too short for the approximation method,
# we would want to not use the approx method



# 2. Epidemic Simulation --------------------------------------------------

## Model parameters
param <- param.net(inf.prob = 0.4, act.rate = 2, rec.rate = 0.01)
# inf.prob - infection probability
# act.rate - given that a partnership exists, how many acts do ppl within
#            that partnership have at the time step
# rec.rate - recovery rate (exponential decay)
# inf.prob + act.rate - deal with the force of infection component


## Initial conditions
init <- init.net(i.num = 10) # 10 people are infected to begin with

## Control settings
control <- control.net(type = "SIS", nsims = 5, nsteps = 500, verbose.int = 0)
# do 5 simulations over 500 time steps
## Run the network simulation
sim <- netsim(est, param, init, control)

## Print the epidemic model output
sim
# si.flow - flow of people from s to i, i.e. incidence of infections over time

?param.net # gives indepthinfo about the parameters

# 3. Model Analysis -------------------------------------------------------

## Plot the epidemic model output
par(mfrow = c(1, 1))
plot(sim)
# line - mean of simulations
# range - interquartile range

## Demonstrates some different plotting options
par(mfrow = c(1, 2))
plot(sim, sim.lines = TRUE, mean.line = FALSE, qnts = FALSE, popfrac = FALSE)
plot(sim, mean.smooth = FALSE, qnts = 1, qnts.smooth = FALSE, popfrac = FALSE)
# can see individual simulations, and the rough stochasticity (if smoothers are
#       set to false)
# help: ?plot.netsim
?plot.netsim

## Use the y argument to pass different model outcomes
par(mfrow = c(1,1))
plot(sim, y = c("si.flow", "is.flow"), qnts = FALSE,
     ylim = c(0, 10), leg = TRUE, main = "Flow Sizes")

## Static network plots
par(mar = c(0,0,0,0), mfrow = c(1, 2))
plot(sim, type = "network", col.status = TRUE, at = 1, sims = 1)
plot(sim, type = "network", col.status = TRUE, at = 500, sims = 1)
# sims = the simulation number. We simulated 5 times, and here are asking for
#       the first sim

## The summary function for time-specific model outcomes
summary(sim, at = 500)
# at = the timestep. we're asking for the final time step here


## The as.data.frame function to extract epi data
df <- as.data.frame(sim)
head(df, 10)
# there are fractions here because we're summarizing across 5 simulations
# NOTE: useful for customization of plots as needed by you


## Default is to calculate means, use the out argument to get other output
df <- as.data.frame(sim, out = "vals", sim = 5)
head(df, 10)

## Extract individual networkDynamic objects from the simulation object
nw1 <- get_network(sim, sim = 1)
nw1


## The transmission matrix
tm1 <- get_transmat(sim, sim = 1)
head(tm1, 10)
# NOTE: Can use this to create a phylogenic tree - for transmission chain


## Save your data!
save(sim, file = "/path/to/folder/day3-tutorial1-sim.rda")

# igraph - doesn't do dynamic network objects, but you can export static
#       parts of it.
# Also use intergraph

