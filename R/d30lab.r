
##
## Network Modeling for Epidemics 2016
## Day 3 | Lab
##


## Load EpiModel
library(EpiModel)


# 1. Model Estimation -----------------------------------------------------

## Initialize the network
num.m1 <- 361 # Males
num.m2 <- 372 # Females

nw <- network.initialize(n = num.m1 + num.m2, bipartite = num.m1,
                         directed = FALSE)
nw

## Fractional degree distributions by mode
## Mode1 = Male; Mode2 = Female
deg.dist.m1 <- c(0.283, 0.568, 0.128, 0.02)
deg.dist.m2 <- c(0.245, 0.742, 0.011, 0.003)

## Check for balancing of edges across modes
check_bip_degdist(num.m1, num.m2,
                  deg.dist.m1, deg.dist.m2)


## -----------------------------------------------------------------------------
## Adding the balanced one provided to us:
deg.dist.m1 <- c(0.328, 0.524, 0.127, 0.021)
deg.dist.m2 <- c(0.20, 0.786, 0.011, 0.003)

## Check for balancing of edges across modes
check_bip_degdist(num.m1, num.m2,
                  deg.dist.m1, deg.dist.m2)



## -----------------------------------------------------------------------------
## Formation formula

formation <- ~edges + b1degree(0:1) + b2degree(0:1)
target.stats <- c(303.601, 118.408, 189.164, 74.400, 292.392)


## Dissolution coefficient calculations
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 38)
coef.diss

## Estimate the model
est <- netest(nw, formation, target.stats, coef.diss)

## Diagnose the model fix
dx <- netdx(est, nsims = 5, nsteps = 500,
            nwstats.formula = ~edges + b1degree(0:3) + b2degree(0:3))
dx

## Compare the model fit statistics back to the "empirical" degree distribution
check_bip_degdist(num.m1, num.m2,
                  deg.dist.m1, deg.dist.m2)

## Plot the diagnostics for formation formula
plot(dx, stats = c("edges", "b1deg1", "b2deg1"))

## Plot dissolution model diagnostics
plot(dx, type = "duration")



# 2. Epidemic Simulation --------------------------------------------------

## Model parameters by mode: .m2 are for mode 2 persons
param <- param.net(inf.prob = 0.2, inf.prob.m2 = 0.2,
                   rec.rate = 0.02, rec.rate.m2 = 0.02)

## Initial conditions
init <- init.net(i.num = 10, i.num.m2 = 10,
                 r.num = 0, r.num.m2 = 0, status.rand = FALSE)
# can set infectious people only in males, etc. can play around
# have to specify the reco ppl incase there already are rec

## Control settings
control <- control.net(type = "SIS", nsims = 5, nsteps = 500)

## Simulate the model
sim <- netsim(est, param, init, control)



# 3. Model Analysis -------------------------------------------------------

## Plot the simulation
plot(sim)

## Plot mode-specific data on prevalence and incidence
par(mfrow = c(1, 2))
plot(sim, y = c("i.num", "i.num.m2"),
     qnts = 0.5, ylim = c(0, 0.4), leg = TRUE)
plot(sim, y = c("si.flow", "si.flow.m2"),
     qnts = 0.5, ylim = c(0, 3), leg = TRUE)

## Static network plot for "most average" simulation
par(mfrow = c(1,1), mar = c(0,0,0,0))
plot(sim, type = "network", col.status = TRUE, at = 50,
     sims = "mean", shp.bip = "square")
# sims = mean - for getting mean values, not an outlier








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

