load("./essentials_allBRCA.RData") # data to work on
load("/Users/Xinyu/Documents/Bayes/code/methylationID.Rdata")
n1 = names(gbIND$SID)
n2 = names(prIND$SID)
n0 = rownames(counts)[which(rowSums(counts) != 0)]
n_int <- intersect(intersect(n1, n2), n0)
#############------------------------------------------------------##############
#############-----------------     expression      ----------------##############
#############------------------------------------------------------##############
#### variables: counts, factors_ls
## parameter 1. n0 --- the genes to be crossed
n0 = rownames(counts)[which(rowSums(counts) != 0)]
## parameter 2. nn --- sampls
nn = c(G1[1:3], G2[1:3])
## parameters 3. multifier of standard deviation 
k = 2
## parameters 4. number of bins
nb <- 25

# normalization: expression = counts / factors
expr <- mapply(`/`, counts[n0, nn], factors_ls[nn])
#counts0 <- counts[n0, nn]
CPM <- sort(unlist(expr))
breaks <- CPM[seq(nb-1) / nb * length(CPM)]
maxCPM <- max(CPM)
max_boundary <- maxCPM+(4*sqrt(maxCPM))+20 
breaks <- unique(breaks)
if (length(breaks) == 1) breaks <- seq(maxPM,max_boundary-maxCPM,(max_boundary-maxCPM)/24) else if (length(breaks) != (nb-1)){
  breaks <- breaks[breaks != 0]
  breaks <- c(breaks, exp(seq(from = log(max(breaks)), to = log(max_boundary), length.out = nb-1-length(breaks)+2))[-c(1, nb-1-length(breaks)+2)] )
}
breaksEXPRESSION <- c(0,breaks,max_boundary)

## result1
breaksEXPRESSION1 <- c(0.000000e+00,4.091192e-02,1.142843e-01,2.902888e-01,6.582777e-01,1.329842e+00,2.416708e+00,3.964829e+00,6.000822e+00,8.506019e+00,1.136917e+01,1.467737e+01,1.836680e+01,2.263194e+01,2.763255e+01,3.359959e+01,4.087030e+01,5.013933e+01,6.273463e+01,8.094219e+01,1.115939e+02,1.896524e+02,5.446752e+02,1.564288e+03,4.492583e+03,1.290254e+04)
## result use all samples
breaksEXPRESSION_all <- c(0.000000e+00,3.141383e-02,9.679162e-02,2.479091e-01,5.692551e-01,1.182731e+00,2.219819e+00,3.773794e+00,5.855501e+00,8.387344e+00,1.132239e+01,1.467913e+01,1.851976e+01,2.296866e+01,2.820315e+01,3.448617e+01,4.221585e+01,5.198129e+01,6.504883e+01,8.401590e+01,1.162970e+02,1.953982e+02,9.114682e+02,4.251700e+03,1.983279e+04,9.251344e+04)


## method 2: same as PINCAGE
## define a function generate a table of CPM and pdf for each sample in Poisson distribution 
## geneCounts - all gene counts for the sample
## fator      - factor for the sample
## k          - multiplier of sd expanding
sample2pdf <- function(geneCounts, factor, k = 2) {
  tempAN <- matrix(ncol=2, nrow=0)
  colnames(tempAN) <- c("cpm","density")
  for(i in 1:length(geneCounts)) {
    lambda <- geneCounts[i]
    X <- seq(round(max(lambda-(k*sqrt(lambda)),0)),round(lambda+(k*sqrt(lambda)))) 
    tempAN <- rbind(tempAN,cbind(X/factor,dpois(X,lambda=lambda))) 
  }
  return(tempAN)
}

library(parallel)
temps <- mcmapply(sample2pdf, counts0, factors_ls[nn])
tempAN <- Reduce(rbind, temps)
tempAN <- as.data.frame(tempAN)
tempAN <- tempAN[order(tempAN$cpm),]
tempAN[,3] <- cumsum(tempAN[,2])
if (max(tempAN[,3]) != 0) tempAN[,3] <- tempAN[,3]/max(tempAN[,3])
breaks <- NULL
noBreaks <- res_expr-1
for (j in 1:noBreaks) { breaks <- c(breaks, tempAN[which(tempAN[,3] >= j*(1/(1+noBreaks))),1][1])}
max <- max(tempAN[nrow(tempAN), 1])
max_boundary <- 20+max+(4*sqrt(max))
breaks <- unique(breaks)
if (length(breaks) == 1) breaks <- seq(max(tempAN[,1]),max_boundary-max(tempAN[,1]),(max_boundary-max(tempAN[,1]))/24) else if (length(breaks) != noBreaks){
  breaks <- breaks[breaks != 0]
  breaks <- c(breaks, exp(seq(from = log(max(breaks)), to = log(max_boundary), length.out = noBreaks-length(breaks)+2))[-c(1, noBreaks-length(breaks)+2)] )
}
breaksEXPRESSION <- c(0,breaks,max_boundary)

# result2 
breaksEXPRESSION2<-c(0.000000e+00,3.039885e-02,9.676293e-02,2.571781e-01,5.805776e-01,1.209537e+00,2.229700e+00,3.712091e+00,5.697890e+00,8.186836e+00,1.102435e+01,1.430479e+01,1.800084e+01,2.226716e+01,2.723737e+01,3.315911e+01,4.044032e+01,4.968692e+01,6.217067e+01,8.032555e+01,1.108342e+02,1.882153e+02,5.419253e+02,1.560357e+03,4.492711e+03,1.293579e+04)

#############------------------------------------------------------##############
#############-----------------     methylation     ----------------##############
#############------------------------------------------------------##############
s1 <- unlist(gbIND$SID)
s2 <- unlist(prIND$SID)
mmatrix_pc0 <- mmatrix_pc[union(s1,s2), nn]
density <- density(mmatrix_pc0,bw=0.14,from=-7,to=7,n=2801,na.rm=TRUE)
density$y <- density$y/sum(density$y)
density$y <- cumsum(density$y)
breaks <- NULL
noBreaks <- res_gb-1
for (j in 1:noBreaks) { breaks <- c (breaks, density$x[which(density$y >= j*(1/(1+noBreaks)))][1])}
breaksMETHYLATION <- sort(c(-7.01,breaks,7.01))

# result3
breaksMETHYLATION <- c(-7.010,-5.825,-5.475,-5.210,-4.965,-4.710,-4.420,-4.050,-3.535,-2.770,-1.780,-0.880,-0.155,0.480,1.070,1.630,2.150,2.620,3.050,3.440,3.800,4.145,4.485,4.845,5.310,7.010)




