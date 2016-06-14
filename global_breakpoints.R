load("./essentials_allBRCA.RData") # data to work on
#############-----------------     expression      ----------------##############
#### variables: counts, factors_ls
## parameter 1. n0 --- the genes to be crossed
n0 = rownames(counts)[which(rowSums(counts) != 0)]
## parameter 2. nn --- sampls
nn = c(G1[1:3], G2[1:3])
## parameters 3. multifier of standard deviation (only used in original expansion)
k = 2
## parameters 4. number of bins
nb <- 25


# normalization: expression = counts / factors
expr <- mapply(`/`, counts[n0, nn], factors_ls[nn])
counts0 <- counts[n0, nn]
CPM <- sort(unlist(expr))
breaksExpression <- CPM[seq(nb-1) / nb * length(CPM)]

tempAN <- matrix(ncol=2, nrow=0)
colnames(tempAN) <- c("cpm","density")

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
tempAN <- Reduced(rbind, temps)
tempAN <- as.data.frame(tempAN)
tempAN <- tempAN[order(tempAN$cpm),]
tempAN[,3] <- cumsum(tempAN[,2])
if (max(tempAN[,3]) != 0) tempAN[,3] <- tempAN[,3]/max(tempAN[,3])
breaks <- NULL
noBreaks <- res_expr-1
for (j in 1:noBreaks) { breaks <- c(breaks, tempAN[which(tempAN[,3] >= j*(1/(1+noBreaks))),1][1])}
max <- max(temp)
if (max > 0) max_boundary <- 20+max+(4*max*max^(-1/2)) else max_boundary <- 20
breaks <- unique(breaks)
if (length(breaks) == 1) breaks <- seq(max(tempAN[,1]),max_boundary-max(tempAN[,1]),(max_boundary-max(tempAN[,1]))/24) else if (length(breaks) != 24){
  breaks <- breaks[-which(breaks %in% 0)]
  breaks <- c(breaks,seq(max(tempAN[,1])+max(breaks),max_boundary-max(tempAN[,1]),(max_boundary-max(tempAN[,1])-max(breaks))/(24-length(breaks))))
}
breaksEXPRESSION <- c(0,breaks,max_boundary)
