
## Initialize the network
num.m1 <- 361 # Males
num.m2 <- 372 # Females

nw <- network.initialize(n = num.m1 + num.m2, bipartite = num.m1,
                         directed = FALSE)
nw
# specify the bipartite argument to create a bipartite network
# No explicit naming - though assume that m1 consists of females

## Fractional degree distributions by mode
## Mode1 = Male; Mode2 = Female
deg.dist.m1 <- c(0.283, 0.568, 0.128, 0.02)
deg.dist.m2 <- c(0.245, 0.742, 0.011, 0.003)

## Check for balancing of edges across modes
check_bip_degdist(num.m1, num.m2,
                  deg.dist.m1, deg.dist.m2)


df = data.frame(deg = 0:3, m1deg = deg.dist.m1/sum(deg.dist.m1), m2deg = deg.dist.m2/sum(deg.dist.m2))
df

bal = check_bip_degdist(num.m1, num.m2,
                  df$m1deg, df$m2deg)



degm1 = 
degm2 = deg.dist.m2 + ((deg.dist.m2/sum(deg.dist.m2)) * .11)
