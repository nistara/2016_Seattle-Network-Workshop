## ----setup, include=FALSE------------------------------------------------
library(EpiModel)
library(knitr)
knitr::opts_chunk$set(cache=TRUE, comment=NA, fig.width=9,
                      fig.height=5.4,
                      cache.path = 'DocumentName_cache/',
                      fig.path='DocumentName_figure/')
op <- par()
par(mfrow=c(1,1),mar=c(3,3,1,1), mgp=c(2,1,0))
np <- par()

## ------------------------------------------------------------------------
set.seed(0)

## ------------------------------------------------------------------------
library(EpiModel)

## ------------------------------------------------------------------------
net1 <- network.initialize(100, directed = FALSE)

## ------------------------------------------------------------------------
coef.diss.1 <- dissolution_coefs(~offset(edges), 90)
coef.diss.1

## ------------------------------------------------------------------------
fit1 <- netest(net1, 
                       formation = ~edges, 
                       target.stats = 20,
                       coef.diss = coef.diss.1
                       )

## ------------------------------------------------------------------------
summary(fit1)

## ------------------------------------------------------------------------
sim1 <- netdx(fit1, nsteps=1000, nsims = 10,
                      keep.tedgelist = TRUE)

## ------------------------------------------------------------------------
sim1

## ------------------------------------------------------------------------
plot(sim1, type='formation')

## ------------------------------------------------------------------------
plot(sim1, type='duration')

## ------------------------------------------------------------------------
plot(sim1, type='dissolution')

## ------------------------------------------------------------------------
sim1$tedgelist[[1]][1:5,]

## ------------------------------------------------------------------------
tel <- sim1$tedgelist[[1]]
hist(tel$duration)
mean(tel$duration[tel$onset < 100])
sum(tel$terminus.censored==TRUE)
plot(tel$onset, tel$terminus)
table(c(tel$head,tel$tail))
hist(table(c(tel$head,tel$tail)))

## ------------------------------------------------------------------------
net2 <- network.initialize(1000, directed = FALSE)
fit2 <- netest(net2, 
   formation=~edges, 
   target.stats=200,
   coef.diss=dissolution_coefs(~offset(edges), 90)
   )
sim2 <- netdx(fit2, nsteps=1000, nsims=10, keep.tedgelist=TRUE)
plot(sim2, type='formation')
plot(sim2, type='duration')
plot(sim2, type='dissolution')

## ------------------------------------------------------------------------
n <- 500
net3 <- network.initialize(n, directed=FALSE)  
net3 %v% 'race' <- c(rep("B",n/2), rep("W",n/2)) #%v% means assign vertex
## attribute called race to net3 graph
net3
table(net3 %v% 'race')


## ------------------------------------------------------------------------
form.formula.3 <- ~edges+nodematch('race')+degree(0)+
    concurrent   ##(no o fnodes with degree 2 or more)
## don' specify degree(1) coz then they'll add up to no of vertices, which is
## a fixed attribute. So adding this will cause colinearity

## NOTE: This following part is very imp, and you need to think through the
## target stats you're ggoing to feed in
target.stats.3 <- c(0.9*n/2, (0.9*n/2)*(5/6), 0.36*n, 0.18*n)

## ------------------------------------------------------------------------
diss.formula.3 <- ~offset(edges)+offset(nodematch("race"))

## ---- eval=FALSE---------------------------------------------------------
## ?dissolution_coefs

## ----results='hide', warning=FALSE---------------------------------------
fit3 <- netest(net3,
             formation= form.formula.3,
             target.stats= target.stats.3,
             coef.diss = dissolution_coefs(
                  ~offset(edges)+offset(nodematch("race")),
                  c(200, 100)
             )
           )

## ---- warning=FALSE------------------------------------------------------
sim3 <- netdx(fit3, nsteps=1000, nsims=10, keep.tedgelist=TRUE)

## ------------------------------------------------------------------------
sim3

## ------------------------------------------------------------------------
plot(sim3, type='formation')

## ---- eval=FALSE---------------------------------------------------------
## plot(sim3, type='duration')

## ------------------------------------------------------------------------
race <- net3 %v% 'race'
tel3 <- sim3$tedgelist[[1]]
mean(tel3$duration[(race[tel3$tail] !=  race[tel3$head]) & tel3$onset < 100])
mean(tel3$duration[(race[tel3$tail] ==  race[tel3$head]) & tel3$onset < 100])

## ------------------------------------------------------------------------
fit4 <- netest(net3,
             formation= form.formula.3,
             target.stats= target.stats.3,
             coef.diss = dissolution_coefs(
                  ~offset(edges)+offset(nodematch("race")),
                  c(20, 10)
             )
)

## ------------------------------------------------------------------------
sim4 <- netdx(fit4, nsteps=1000, nsims=10, keep.tedgelist=TRUE,)

## ------------------------------------------------------------------------
plot(sim4, type='formation')

## ------------------------------------------------------------------------
fit5 <- netest(net3,
             formation= form.formula.3,
             target.stats= target.stats.3,
             coef.diss = dissolution_coefs(
                  ~offset(edges)+offset(nodematch("race")),
                  c(20, 10)
             ),
             edapprox=FALSE
           )
## edapprox - edge duration approximation
## Now you're running a STRGEM for the first time actually!

## ---- warning=FALSE, message=FALSE---------------------------------------
sim5 <- netdx(fit5, nsteps=1000, nsims=10, keep.tedgelist=TRUE)

## ------------------------------------------------------------------------
plot(sim5, type='formation')

## ------------------------------------------------------------------------
race <- net3 %v% 'race'
tel5 <- sim5$tedgelist[[1]]
mean(tel5$duration[(race[tel5$tail] !=  race[tel5$head]) & tel5$onset<100])
mean(tel5$duration[(race[tel5$tail] ==  race[tel5$head]) & tel5$onset<100])

