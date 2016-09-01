
##
## Network Modeling for Epidemics 2015
## Day 4 | Tutorial 1: Dynamic Network Visualization
##

library("EpiModel")
library("ndtv")

## Basic model estimation and simulation
nw <- network.initialize(n = 100, directed = FALSE)
formation <- ~edges + concurrent
target.stats <- c(45, 25)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 5)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
param <- param.net(inf.prob = 1)
init <- init.net(i.num = 1)
control <- control.net(type = "SI", nsteps = 25, nsims = 1, verbose = FALSE)
sim <- netsim(est, param, init, control)

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
