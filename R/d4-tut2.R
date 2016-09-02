
##
## Network Modeling for Epidemics 2015
## Day 4 | Tutorial 2: Demography and Attributes
##

library(EpiModel)

## Initialize the network
nw <- network.initialize(n = 500, directed = FALSE)
nw <- set.vertex.attribute(nw, attrname = "race", value = rep(0:1, each = 250))

## %v% is the same as set.vertex.attribute
## bipartite | assortative | dissortative mixing

## Formation formulas
formation <- ~edges + nodefactor("race") + nodematch("race")
## nodefactor and match = attribute based
## nfactor - controlling differential mean degree by attri
## match- level of between grp mixing?


## Target stats for model 1
target.stats <- c(125, 125, 62.5)

## Dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 40, d.rate = 0.01)
coef.diss
## i.e. death rate = 1% per time step
## Coefficient - represents the intensity of persistence (adjusted is higher than crude)

## Estimate model 1
est1 <- netest(nw, formation, target.stats, coef.diss)

## Summary of model fit
summary(est1)

## Diagnostics on Model 1 with multi-core processing
dx1 <- netdx(est1, nsims = 10, nsteps = 1000, ncores = 4,
             nwstats.formula = ~edges +
                                nodefactor("race", base = 0) +
                                nodematch("race"))
dx1
## ncores allows you to run the simulations in parallel. It could take over you computer
## if you use all coresin your laptop. But it increases efficiency
## base = 0 - means don't use any one as a reference


plot(dx1)

## Target stats for model 2
target.stats <- c(125, 187.5, 112.5)

## Estimate model 2
est2 <- netest(nw, formation, target.stats, coef.diss)
summary(est2)

## Diagnostics for Model 2
dx2 <- netdx(est2, nsims = 10, nsteps = 1000, ncores = 4,
            nwstats.formula = ~edges +
                               nodefactor("race", base = 0) +
                               nodematch("race"))
dx2
## if you remove base = 0, it treats 0 as the reference. HEre you want to see the target stats for all of them. Hence we include base=0. Try without it.
plot(dx2)


## Epidemic model parameters
param <- param.net(inf.prob = 0.1, act.rate = 5,
                   b.rate = 0.01, ds.rate = 0.01, di.rate = 0.01)
## this death rate should be the same as one specified in the dissolution - .01
## if mortality is diff - take weighted average (calibrated across time as well)
## Keep females the first node - coz birth rate is a function of females

## 1% f sus and inf die at every tie step
init <- init.net(i.num = 50)
control <- control.net(type = "SI", nsteps = 500, nsims = 1, epi.by = "race")
## Just doing one simulation for now. oirignal was 10

## Simulate the two models
sim1 <- netsim(est1, param, init, control)
sim2 <- netsim(est2, param, init, control)

## Show model output
sim1

## Plot the edges fit after epidemic simulation
par(mfrow = c(1, 2))
plot(sim1, type = "formation", stats = "edges")
plot(sim2, type = "formation", stats = "edges")

dev.off()
## Plot the overall prevalence from each model
plot(sim1, y = "i.num", qnts = 1, mean.col = "steelblue",
     qnts.col = "steelblue", main = "Total Prevalence")
plot(sim2, y = "i.num", qnts = 1, mean.col = "firebrick",
     qnts.col = "firebrick", add = TRUE)
legend("topleft", c("Model 1", "Model 2"), lwd = 3,
       col = c("steelblue", "firebrick"), bty = "n")
## add = TRUE - adds the second plot to the first one


## Plot the race-specific prevalence from each model
par(mfrow = c(1, 2))
plot(sim1, y = c("i.num.race0", "i.num.race1"),  leg = TRUE, qnts = 1,
     ylim = c(0, 200), main = "M1: Disease Prevalence by Race")
plot(sim2, y = c("i.num.race0", "i.num.race1"), leg = TRUE,  qnts = 1,
     ylim = c(0, 200), main = "M2: Disease Prevalence by Race")

## Save the data
save(sim1, sim2, file = "path/to/folder/d4-tutorial1.rda")
