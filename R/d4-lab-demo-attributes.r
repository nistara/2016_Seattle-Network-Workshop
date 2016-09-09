## NOTE: because we 


##
## Network Modeling for Epidemics 2016
## Day 4 | Lab
##


## Load EpiModel
library(EpiModel)


# 1. Model Estimation -----------------------------------------------------

## Flip the first mode to females (coz birth rate is calculated based on first mode,
## and ingeneral, it's cal based on females
## Initialize the network
num.m1 <- 372 # Females 
num.m2 <- 361 # Males
totpop = num.m1 + num.m2

nw <- network.initialize(n = num.m1 + num.m2, bipartite = num.m1,
                         directed = FALSE)
## race_value = sample(totpop, num.m1)
## race_num = ifelse(1:totpop %in% race_value, 0, 1)

race_type = sample(0:1, totpop, replace = TRUE)
nw <- set.vertex.attribute(nw, attrname = "race", value = race_type)
nw

## Fractional degree distributions by mode
## Mode1 = Male; Mode2 = Female
deg.dist.m1 <- c(0.155, 0.631, 0.111, 0.103)
deg.dist.m2 <- c(0.173, 0.578, 0.128, 0.121)

## Check for balancing of edges across modes
check_bip_degdist(num.m1, num.m2,
                  deg.dist.m1, deg.dist.m2)




## -----------------------------------------------------------------------------
## Formation formula
## meandeg = edges/(tot_pop/2)
md = 432 /(totpop/2)
race1 = sum(race_type)

## Random mixing with same degree-----------------------------------------------
formation <- ~edges + nodefactor("race") + nodematch("race")
target.stats <- c(432, md*race1, 432/2)


## Dissolution coefficient calculations
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 38,
                               d.rate = 0.01)
coef.diss

## Estimate the model
est <- netest(nw, formation, target.stats, coef.diss)
summary(est)

## Diagnose the model fix
dx <- netdx(est, nsims = 5, nsteps = 100, ncores = 4,
            nwstats.formula = ~edges + nodefactor("race", base = 0) + nodematch("race"))
dx


## Plot the diagnostics for formation formula
plot(dx)
plot(dx, stats = c("edges", "b1deg0", "b2deg0"))

## Plot dissolution model diagnostics
plot(dx, type = "duration")

################################################################################
## Random mixing with different degree------------------------------------------
## increase mean degree 1.5 times
md2 = md*1.5
formation2 <- ~edges + nodefactor("race") + nodematch("race")
target.stats2 <- c(432, md2*race1, 432/2)


## Dissolution coefficient calculations
coef.diss2 <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 38,
                               d.rate = 0.01)
coef.diss2

## Estimate the model
est2 <- netest(nw, formation2, target.stats2, coef.diss2)
summary(est2)
## Diagnose the model fix
dx2 <- netdx(est2, nsims = 5, nsteps = 100, ncores = 4,
            nwstats.formula = ~edges + nodefactor("race", base = 0) + nodematch("race"))
dx2


## Plot the diagnostics for formation formula
plot(dx2)
################################################################################

## Same degree with non-random mixing------------------------------------------
# 80% of the partnerships occur within the same race
formation3 <- ~edges + nodefactor("race") + nodematch("race")
target.stats3 <- c(432, md*race1, 432*.80)


## Dissolution coefficient calculations
coef.diss3 <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 38,
                               d.rate = 0.01)
coef.diss3

## Estimate the model
est3 <- netest(nw, formation3, target.stats3, coef.diss3)
summary(est3)
## Diagnose the model fix
dx3 <- netdx(est3, nsims = 5, nsteps = 100, ncores = 4,
            nwstats.formula = ~edges + nodefactor("race", base = 0) + nodematch("race"))
dx3


## Plot the diagnostics for formation formula
plot(dx3)
################################################################################
## Different degree with non-random mixing------------------------------------------
# 80% of the partnerships occur within the same race
formation4 <- ~edges + nodefactor("race") + nodematch("race")
target.stats4 <- c(432, md*1.5*race1, 432*.80)


## Dissolution coefficient calculations
coef.diss4 <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 38,
                               d.rate = 0.01)
coef.diss4

## Estimate the model
est4 <- netest(nw, formation4, target.stats4, coef.diss4)
summary(est4)
## Diagnose the model fix
dx4 <- netdx(est4, nsims = 5, nsteps = 100, ncores = 4,
            nwstats.formula = ~edges + nodefactor("race", base = 0) + nodematch("race"))
dx4


## Plot the diagnostics for formation formula
plot(dx4)
################################################################################

# 2. Epidemic Simulation --------------------------------------------------

## Model parameters by mode: .m2 are for mode 2 persons
param <- param.net(inf.prob = 0.2, inf.prob.m2 = 0.2,
                   rec.rate = 0.02, rec.rate.m2 = 0.02,
                   b.rate = 0.005, b.rate.m2 = NA,
                   ds.rate = 0.001, ds.rate.m2 = 0.001,
                   di.rate = 0.001, di.rate.m2 = 0.001)

## Initial conditions
init <- init.net(i.num = 10, i.num.m2 = 10,
                 r.num = 0, r.num.m2 = 0, status.rand = FALSE)
# can set infectious people only in males, etc. can play around
# have to specify the reco ppl incase there already are rec

## Control settings
control <- control.net(type = "SIS", nsims = 25, nsteps = 250, epi.by = "race")

## Simulate the model
sim <- netsim(est, param, init, control)
sim1 <- netsim(est2, param, init, control)
sim2 <- netsim(est3, param, init, control)
sim3 <- netsim(est4, param, init, control)

## Note the discrepancy in the names of the sims and the ests. Am leaving as is
## for now

# 3. Model Analysis -------------------------------------------------------

## Plot the simulation
plot(sim)
plot(sim1)
plot(sim2)
plot(sim3)




################################################################################
## Plot the edges fit after epidemic simulation
par(mfrow = c(2, 2))
plot(sim, type = "formation", stats = "edges")
plot(sim1, type = "formation", stats = "edges")
plot(sim2, type = "formation", stats = "edges")
plot(sim3, type = "formation", stats = "edges")
dev.off()


## Plot the overall prevalence from each model
plot(sim, y = "i.num", qnts = 1, mean.col = "steelblue",
     qnts.col = "steelblue", main = "Total Prevalence")
plot(sim1, y = "i.num", qnts = 1, mean.col = "firebrick",
     qnts.col = "firebrick", add = TRUE)
plot(sim2, y = "i.num", qnts = 1, mean.col = "green",
     qnts.col = "firebrick", add = TRUE)
plot(sim3, y = "i.num", qnts = 1, mean.col = "purple",
     qnts.col = "firebrick", add = TRUE)
legend("topleft", c("Model 0", "Model 1", "Model 2", "Model 3"), lwd = 3,
       col = c("steelblue", "firebrick", "green", "purple"), bty = "n")
## add = TRUE - adds the second plot to the first one


## Plot the race-specific prevalence from each model
par(mfrow = c(2, 2))
plot(sim, y = c("i.num.race0", "i.num.race1"),  leg = FALSE, qnts = 1,
     ylim = c(0, 500), main = "M0:Random mixing with same mean degree")
plot(sim1, y = c("i.num.race0", "i.num.race1"), leg = TRUE,  qnts = 1,
     ylim = c(0, 500), main = "M1: Random mixing + higher mean degree")
plot(sim2, y = c("i.num.race0", "i.num.race1"), leg = TRUE,  qnts = 1,
     ylim = c(0, 500), main = "M2: Assortative mixing + same mean degree")
plot(sim3, y = c("i.num.race0", "i.num.race1"), leg = TRUE,  qnts = 1,
     ylim = c(0, 500), main = "M3: ssortative mixing + differential degree")





################################################################################








## as.data.frame function for data extraction to calculate cumulative incidence
df <- as.data.frame(sim)
sum(df$si.flow)
sum(df$si.flow.m2)

## How to hand-calculate the incidence rate per 100 person-years at risk
#Assuming timestep is 1 week
ir.m1 <- (df$si.flow/df$s.num) * 100 * 52
which.max(ir.m1)
ir.m1[which.max(ir.m1)]

## Extract the transmission matrix
tm <- get_transmat(sim)
head(tm, 15)

# at - timestep at which transm occured
# 
## Use the transmission matrix to calculate number of "early stage infections"
mean(tm$infDur < 10)
# "across" time step - hence the infDur

## Save the Data!!
save(sim, file = "/path/to/folder/day3-tutorial2-sim.rda")




## Network Movie----------------------------------------------------------------
## Basic model estimation and simulation
## nw <- network.initialize(n = 100, directed = FALSE)
## formation <- ~edges + concurrent
## target.stats <- c(45, 25)
## coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 5)
## est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
## param <- param.net(inf.prob = 1)
## init <- init.net(i.num = 1)
## control <- control.net(type = "SI", nsteps = 25, nsims = 1, verbose = FALSE)
## sim <- netsim(est, param, init, control)

library(ndtv)

## Extract the network
nw <- get_network(sim)

## Add a time-varying color attribute to vertices
nw <- color_tea(nw, verbose = FALSE)

## Animation options
slice.par <- list(start = 1, end = 25, interval = 1,
                  aggregate.dur = 1, rule = "any")
render.par <- list(tween.frames = 10, show.time = FALSE)
plot.par <- list(mar = c(0, 0, 0, 0))

## Compute the animation positions of nodes
compute.animation(nw, slice.par = slice.par)

## Render the movie
render.d3movie(
    nw,
    render.par = render.par,
    plot.par = plot.par,
    vertex.cex = 0.9,
    vertex.col = "ndtvcol",
    edge.col = "darkgrey",
    vertex.border = "lightgrey",
    displaylabels = FALSE,
    filename = paste0(getwd(), "/movie.html"))

proximity.timeline(nw, vertex.col = "ndtvcol", spline.style = "color.attribute")

plot(sim, type = "network", col.status = TRUE, at = 5)







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
