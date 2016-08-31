
##
## Network Modeling for Epidemics 2016
## Day 3 | Tutorial 4: Interventions
##

## Load EpiModel
library(EpiModel)


## Fit a basic temporal network model
nw <- network.initialize(n = 100, directed = FALSE)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

## Simulate an SI epidemic with highly effective vaccine intervention
param <- param.net(inf.prob = 0.5, inter.eff = 0.96, inter.start = 25)
init <- init.net(i.num = 5)
control <- control.net(type = "SI", nsteps = 100, nsims = 10)
sim <- netsim(est, param, init, control)
plot(sim)

## Simulate an SIS epidemic with condom provision intervention
param <- param.net(inf.prob = 0.5, inter.eff = 0.8, inter.start = 100,
                   rec.rate = 0.07)
init <- init.net(i.num = 10)
control <- control.net(type = "SIS", nsteps = 250, nsims = 10)
sim <- netsim(est, param, init, control)
plot(sim)
