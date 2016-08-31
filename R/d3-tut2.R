
##
## Network Modeling for Epidemics 2016
## Day 3 | Tutorial 2: SIR Epidemic in a Bipartite Network
##

## -----------------------------------------------------------------------------
## NOTE: DCMs and balance - can only balance at the outset, but if the
#       no of males and females changes across time, the model won't
#       be balanced...and we have a silent error.
#       This can be balanced in the network model - pg 8 d3-s1.pdf (with the
#       current sex ratio, and the dominant contact rate - e.g. base male contact
#       rate on that of female contact rate and the sex ratio)
#
# Balance issue - the total number of ties has to be the same (i.e. the no of part-
#       nerships. The mean degree doesn't have to be the same
## -----------------------------------------------------------------------------

## Load EpiModel
library(EpiModel)


# 1. Model Estimation -----------------------------------------------------

## Initialize the network
num.m1 <- num.m2 <- 250
nw <- network.initialize(n = num.m1 + num.m2, bipartite = num.m1,
                         directed = FALSE)
nw
# specify the bipartite argument to create a bipartite network
# No explicit naming - though assume that m1 consists of females

## Fractional degree distributions by mode
## Mode1 = Female; Mode2 = Male
deg.dist.m1 <- c(0.40, 0.55, 0.04, 0.01)
deg.dist.m2 <- c(0.54, 0.31, 0.10, 0.05)
# They add up to 100%
# for 0, 1, 2, and 3 partners
# Assume that mean degree is .66 for both men and women to see
#       non sex-specific degree distribution
# The higher degree men are balanced out by the lower degree men as well

## The expectations under a Poisson distribution
pois.dists <- c(dpois(0:2, lambda = 0.66),
                ppois(2, lambda = 0.66, lower.tail = FALSE))

## Plot the different degree distributions
par(mar = c(3,3,2,1), mgp = c(2,1,0), mfrow = c(1,1))
barplot(cbind(deg.dist.m1, deg.dist.m2, pois.dists),
        beside = TRUE, ylim = c(0, 0.6), col = rainbow(4))
legend("topright", legend = paste0("deg", 3:0),
       pch = 15, col = rev(rainbow(4)),
       cex = 0.9, bg = "white")

## Check for balancing of edges across modes
check_bip_degdist(num.m1, num.m2,
                  deg.dist.m1, deg.dist.m2)
# check bipartite degree distribution
#1st ad 3rd col - fractional deg distr
# 2nd and 4th - counts
# final line - should be balanced - if you're doing heterosexual model

## *
## -----------------------------------------------------------------------------
## Formation formula
formation <- ~edges + b1degree(0:1) + b2degree(0:1)
target.stats <- c(165, 100, 137.5, 135, 77.5)
# allow mode specific statistic to be entered in
# 165 from total edges in table above
# b1 = mode 1
# b2 - mode 2
# b1 degree 0 to 1
# Not passing in degress 2 and 3 in this model right now -explained later
## -----------------------------------------------------------------------------

## Dissolution coefficient calculations
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
coef.diss
# why do we pass the offset term???

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
# shows the burn-in

## *
## WHY DO WE LEAVE OUT 4 parts? the degress 2 and 3?
# This depends on the constraints of number of males/females and their edges
# Chk out picture



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
control <- control.net(type = "SIR", nsims = 5, nsteps = 500)

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

