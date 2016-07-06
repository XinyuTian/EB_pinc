## small sample with global PINC
library(aws)
load("~/Downloads/essentials_allBRCA.RData")
load("/Users/Xinyu/Documents/Bayes/code/methylationID.Rdata") # gbIND prIND
counts0 = counts
mmatrix_pc0 = mmatrix_pc
G10 <- G1
G20 <- G2
## parameter 1. nn --- sampls
G1 <- G10[1:5]
G2 <- G20[1:5]
nn = c(G1, G2)
counts = counts0[, nn]
## parameter 2. n0 --- the genes to be crossed
n0 = rownames(counts)[apply(counts, 1, function(xx) sum(xx != 0) > 1)]
counts = counts[n0, ]
## parameters 3. number of bins
res_expr <- 25
epsilon <- 1e-6

## ------- working list of BRCA --------
n1 = names(gbIND$SID)
n2 = names(prIND$SID)

########################################################################
############ ---------     EXPRESSION    ----------- ###################
########################################################################
# normalization: expression = counts / factors
expr <- mapply(`/`, counts, factors_ls[nn]); rownames(expr) <- n0

# borrow information from other genes in each sample
counts_int <- round(counts)
t1 = counts_int[,1]
t2 = table(t1)
G = length(n0)
h0 <- 1
NN <- function(xx) eval(parse(text = paste0('t2["', xx, '"]')))
pz <- function(z, h=h0) {
  res <- sapply(seq(min(170,z)), function(j) NN(z-j+1) * h^(j-1) * exp(h) / factorial(j-1))
  return(sum(res) / G)
}
delta1 <- function(z, h=h0) {
  return((z+1) * pz(z+1, h=h) / pz(z, h=h) - h)
}
delta2 <- function(y, h=h0) {
  res <- 0
  j <- 1
  addj <- 1
  while(addj > 1e-6) {
    addj <- h^j * exp(-h) / factorial(j) * delta1(y+j, h=h)
    res <- res + addj
    j <- j + 1
  }
  return(res)
}

## global breakpoints for EXPRESSION
CPM <- sort(unlist(normExpr))
breaks <- CPM[seq(res_expr-1) / res_expr * length(CPM)]
max <- max(normExpr)
max_boundary <- 20+max+(4*sqrt(max))
breaks <- unique(breaks)
breaks <- breaks[breaks != 0]
breaks <- c(breaks, exp(seq(from = log(max(breaks)), to = log(max_boundary), length.out = res_expr-length(breaks)+1))[-1] )
breaks <- c(0,breaks)
#  breaks <- c(0.00000000,0.04359852,0.15889635,0.27820479,0.39144912,0.49485813,0.59039591,0.68232976,0.76975281,0.85710496,0.94471197,1.03308237,1.12337415,1.21859973,1.31688533,1.42129337,1.53532082,1.65981629,1.80272405,1.96986175,2.17580121,2.44279177,2.83111287,3.48316340,12.42933771,44.35291085)

## discretization: integration in expression
integrand_e <- function(x,k) {dpois(k,x)}
integrationExpr <- array(dim = c(length(n0), length(nn), res_expr), dimnames = list(n0, nn, 1:res_expr))
for (current_gene in n0) {
  current_gene_breaks <- breaks * meanExpr[current_gene]
  for (current_sample in nn) {
    current_breaks <- current_gene_breaks * factors_ls[current_sample]
    read_count <- round(counts[current_gene, current_sample])
    for (freq in 1:res_expr) {
      integrationExpr[current_gene, current_sample, freq] <- integrate(integrand_e, lower = current_breaks[freq], upper = current_breaks[freq+1], read_count,stop.on.error=FALSE)$value
    }
  }
}


########################################################################
############ ---------     METHYLATION    ----------- ##################
########################################################################
# s1 <- unlist(gbIND$SID)
# s2 <- unlist(prIND$SID)
# mmatrix_pc <- mmatrix_pc0[union(s1,s2), nn]
# CpG <- sort(c(mmatrix_pc))
# breaks <- CpG[seq(res_expr-1) / res_expr * length(CpG)]
# breaksMETHYLATION <- sort(c(-7.01,breaks,7.01))
# # breaksMETHYLATION0 <- c(-7.0100000,-5.8478940,-5.4713917,-5.1926302,-4.9391243,-4.6778994,-4.3808218,-4.0109015,-3.4968734,-2.7292534,-1.7247853,-0.8172996,-0.1076732,0.5003966,1.0753891,1.6241419,2.1338372,2.6021923,3.0282594,3.4149241,3.7716902,4.1075183,4.4396423,4.8021497,5.2720302,7.0100000)
# 
n.meth <- intersect(n1, n2)
s1 <- unlist(gbIND$SID[n.meth])
s2 <- unlist(prIND$SID[n.meth])
ns <- union(s1,s2)
mmatrix_pc <- mmatrix_pc0[ns, nn]
CpG <- sort(c(mmatrix_pc))
breaks <- CpG[seq(res_expr-1) / res_expr * length(CpG)]
breaksMETHYLATION <- sort(c(-7.01,breaks,7.01))
# breaksMETHYLATION <- c(-7.0100000,-5.8518614,-5.4764971,-5.1987799,-4.9469305,-4.6881684,-4.3951446,-4.0319128,-3.5305775,-2.7834195,-1.7849408,-0.8625267,-0.1391094,0.4759792,1.0564102,1.6100164,2.1235138,2.5953808,3.0243073,3.4126582,3.7709041,4.1075320,4.4404927,4.8034401,5.2734790,7.0100000)

## discretization: integration in methylation
## Gaussian integrand, sd = 0.14
integrand_m <- function(x,mean) {dnorm(x=mean,mean=x,sd=0.14)}
integrationMeth <- array(dim = c(length(ns), length(nn), res_expr), dimnames = list(ns, nn, 1:res_expr))
for (current_site in ns) {
  for (current_sample in nn) {
    current_CpG <- mmatrix_pc[current_site,current_sample]
    for (freq in 1:res_expr) {
      integrationMeth[current_site, current_sample, freq] <- integrate(integrand_m, lower = breaksMETHYLATION[freq], upper = breaksMETHYLATION[freq+1], current_CpG,stop.on.error=FALSE)$value
    }
  }
}
for (current_site in ns) {
  i=i+1; if(i%%5000 == 0) print(i)
  for (current_sample in nn) {
    integrationMeth[current_site, current_sample, ] <- (integrationMeth[current_site, current_sample, ] + epsilon) / sum(integrationMeth[current_site, current_sample, ] + epsilon)
  }
}

