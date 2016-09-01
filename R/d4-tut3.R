
##
## Network Modeling for Epidemics 2015
## Day 4 | Tutorial 3: Serosorting
##

library(EpiModel)


## Initialize network
n <- 500
nw <- network.initialize(n, directed = FALSE)

## Set status attribute on network
prev <- 0.2
infIds <- sample(1:n, n*prev)
nw <- set.vertex.attribute(nw, "status", "s")
nw <- set.vertex.attribute(nw, "status", "i", infIds)

## Degree by infection status
mean.deg.inf <- 0.3
inedges.inf <- mean.deg.inf * n * prev

mean.deg.sus <- 0.8
inedges.sus <- mean.deg.sus * n * (1 - prev)

## Total edges
edges <- (inedges.inf + inedges.sus)/2

## Hardy Weinberg method for conditional probability
p <- inedges.sus/(edges*2)
q <- 1 - p
nn <- p^2
np <- 2*p*q
pp <- q^2
round(nn + pp, 3)

## Simulation method for probabilities
fit <- ergm(nw ~ edges + nodefactor("status"),
            target.stats = c(edges, inedges.sus))
sim <- simulate(fit, statsonly = TRUE, nsim = 1e5,
                monitor = ~nodematch("status"))
mn <- colMeans(sim)
mn
round(as.numeric(mn[3]/mn[1]), 3)

## Our nodematch statistic
nmatch <- edges * 0.91

## Model 1 estimation
formation <- ~edges + nodefactor("status") + nodematch("status")
target.stats <- c(edges, inedges.sus, nmatch)

coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 50)

est <- netest(nw, formation, target.stats, coef.diss)

## Model 1 diagnostics
dx <- netdx(est, nsims = 5, nsteps = 500,
            nwstats.formula = ~edges +
                               meandeg +
                               nodefactor("status", base = 0) +
                               nodematch("status"), verbose = FALSE)
dx
plot(dx, plots.joined = FALSE)
plot(dx, type = "duration")


## Model 2 estimation and diagnostics
nw <- network.initialize(n, directed = FALSE)
est2 <- netest(nw, formation = ~edges, target.stats = edges, coef.diss)

dx2 <- netdx(est2, nsims = 5, nsteps = 1000,
             nwstats.formula = ~edges + meandeg)
dx2
plot(dx2, plots.joined = FALSE)

## Model 1 epidemic simulation
param <- param.net(inf.prob = 0.03)
init <- init.net()
control <- control.net(type = "SI", nsteps = 500, nsims = 5,
                       nwstats.formula = ~edges +
                                          meandeg +
                                          nodefactor("status", base = 0) +
                                          nodematch("status"),
                       save.network = FALSE)

sim <- netsim(est, param, init, control)


## Model 2 epidemic simulation
param <- param.net(inf.prob = 0.03)
init <- init.net(i.num = n*prev, status.rand = FALSE)
control <- control.net(type = "SI", nsteps = 500, nsims = 5,
                       nwstats.formula = ~edges + meandeg,
                       save.network = FALSE)
sim2 <- netsim(est2, param, init, control)


## Compare epidemics
par(mfrow = c(1,2))
plot(sim, main = "Serosorting")
plot(sim2, main = "No Serosorting")

## Overlaid plot
par(mfrow = c(1,1))
plot(sim, y = "i.num", sim.lines = FALSE, qnts = 1)
plot(sim2, y = "i.num", sim.lines = FALSE, qnts = 1,
     mean.col = "firebrick", qnts.col = "firebrick", add = TRUE)
legend("topleft", c("Serosort", "Non-Serosort"), lty = 1, lwd = 3,
       col = c("steelblue", "firebrick"), cex = 0.9, bty = "n")

## Network statistics
plot(sim, type = "formation", plots.joined = FALSE)

plot(sim, type = "formation", stats = c("nodefactor.status.s",
                                        "nodefactor.status.i"))

plot(sim, type = "formation",
     stats = c("edges", "nodematch.status"),
     sims = 1, sim.lwd = 2)
