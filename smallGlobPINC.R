## small sample with global PINC
library(aws)
load("essentials_allBRCA.RData")
counts0 = counts
## parameter 1. nn --- sampls
nn = c(G1[1:5], G2[1:5])
counts = counts0[, nn]
## parameter 2. n0 --- the genes to be crossed
n0 = rownames(counts)[apply(counts, 1, function(xx) sum(xx != 0) > 1)]
counts = counts[n0, ]
## parameters 3. number of bins
res_expr <- 25

mmatrix_pc0 = mmatrix_pc
mmatrix_pc = mmatrix_pc[, nn]
G10 <- G1
G20 <- G2
G1 <- G10[1:5]
G2 <- G20[1:5]

# normalization: expression = counts / factors
expr <- mapply(`/`, counts, factors_ls[nn]); rownames(expr) <- n0
## ------- standardization --------
## do not borrow information from means, only borrow information from sd
## as a result, scale mean = 0, sd = shrinked (sd / mean) 
## shrinkage = smooth: kernsm

# function for normalization
normCoefSD <- function(values, h=1) {
  sds <- apply(values, 1, sd)
  if(any(sds == 0)) {
    rm_ind <- which(sds == 0)
    sds <- sds[-rm_ind]
    values <- values[-rm_ind, ]
  }
  means <- apply(values, 1, mean)
  prior_valuesSD <- kernsm(sds / means, h=h)
  sdsSM <- prior_valuesSD@yhat
  values <- (values - means) / sds * c(sdsSM)
  return(list("values" = values, 'sdSM' = sdsSM))
}

normExprRes <- normCoefSD(expr)
normExpr <- normExprRes$values
sdExprSM <- normExprRes$sdSM

exprGroup <- rbind(expr[, 1:5], expr[, 6:10])

## global breakpoints for EXPRESSION
CPM <- sort(unlist(normExpr))
breaks <- CPM[seq(res_expr-1) / res_expr * length(CPM)]
breaksEXPRESSION <- c(-7,breaks,-7)
## discretization: integration in expression
frequencies_expr <- rep(0,length(breaksEXPRESSION)-1)
for (current_gene in n0) {
  sd_current_gene <- sdExprGenesSM[current_gene]
  for (current_sample in nn) {
    read_expr <- expr[current_gene, current_sample]
    
  }
}
for (freq in 1:res_expr) {
  frequencies_expr[freq] <- integrate(integrand_m, lower = breaksEXPRESSION[freq], upper = breaksEXPRESSION[freq+1], expr,stop.on.error=FALSE)[1]
}
frequencies_expr <- unlist(frequencies_expr)
if (all(frequencies_expr == 0)) frequencies_expr[length(frequencies_expr)] <- 1
frequencies_expr <- frequencies_expr + epsilon_e_G1
frequencies_expr <- frequencies_expr/sum(frequencies_expr)
