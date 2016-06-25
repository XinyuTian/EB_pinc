## GOES WITH calLike FUNCTION

load("/Users/Xinyu/Documents/Bayes/code/integrationExpr.Rdata") 
load("/Users/Xinyu/Documents/Bayes/code/integrationMeth.Rdata") 
load("/Users/Xinyu/Documents/Bayes/code/methylationID.Rdata") # gbIND prIND
load("/Users/Xinyu/Documents/Bayes/code/workingList.Rdata") # workingList, G1, G2, template_promoter_CpGs, template_body_CpGs, nruns
library(smoothie)
library(aws)

res_expr <- dim(integrationExpr)[3]

tensor_product <- function(matrix1,matrix2,smooth_h=0,normalize=c("row","column","no"),kernel=c("gauss","cauchy","minvar")) {
  if (is.matrix(matrix1) && is.matrix(matrix2)) result <- matrix(ncol=ncol(matrix1),nrow=ncol(matrix2),data=rep(0,ncol(matrix1)*ncol(matrix2))) else result <- matrix(ncol=length(matrix1),nrow=length(matrix1),data=rep(0,length(matrix1)*length(matrix2)))
  
  if (is.matrix(matrix1) && is.matrix(matrix2)) for (i in 1:nrow(matrix1)) {
    result <- result + matrix(ncol=ncol(matrix1),nrow=ncol(matrix2),byrow=TRUE,data=apply(expand.grid(matrix1[i,],matrix2[i,]), 1, prod))
  } else result <- result + matrix(nrow=length(matrix1),ncol=length(matrix2),byrow=TRUE,data=apply(expand.grid(matrix1,matrix2), 1, prod))
  
  if (is.matrix(matrix1) && is.matrix(matrix2)) result <- result/nrow(matrix1)
  if (!is.null(kernel) && smooth_h > 0) result <- kernel2dsmooth(result,kernel.type = kernel[1],sigma=smooth_h,nx=ncol(matrix2),ny=ncol(matrix1))
  if (normalize[1] == "row") for (i in 1:nrow(result)) result[i,] <- result[i,]/sum(result[i,]) else if (normalize[1] == "column") for (i in 1:nrow(result)) result[,i] <- result[,i]/sum(result[,i])
  return(result)
}
geo_mean <- function(data) {
  log_data <- log(data)
  gm <- exp(mean(log_data[is.finite(log_data)]))
  return(gm)
}

# function to generate potential file
pot <- function(prGeneCpG, gbGeneCpG, exprGene, sm_param) {
  smooth_e <- sm_param[1]
  smooth_gb <- sm_param[2]
  smooth_pr <- sm_param[3]
  
  prior_expr <- kernsm(apply(exprGene,2,mean),h=smooth_e)
  prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
  
  M_exprGB <- tensor_product(gbGeneCpG,exprGene,smooth_h=smooth_gb)
  M_exprP <- tensor_product(prGeneCpG,exprGene,smooth_h=smooth_pr)
  return(list = list(M_exprGB, M_exprP, prior_expr))
}

# function to calculate likelihood
calLike <- function(exprGene, integrationMeth, prior) {
  result = rep(0, nrow(exprGene))
  for(i in 1:nrow(exprGene)) {
    EXPR.likelihood = exprGene[i, ] 
    mP.likelihood = integrationMeth[IDs_promoter, i, , drop=F]
    mGB.likelihood = integrationMeth[IDs_body, i, , drop=F]
    
    M_exprGB <- prior[[1]]
    M_exprP <- prior[[2]]
    expr <- prior[[3]]
    sample1_G1model_liks = 0
    for (j in seq(25)){
      sum_k = 0 
      for (k in seq(25)) {
        sum_k = sum_k + M_exprP[j,k] * prod(mP.likelihood[, 1, k])
      }
      sum_l = 0
      for (l in seq(25)) {
        sum_l = sum_l + M_exprGB[j,l] * prod(mGB.likelihood[, 1, l])
      }
      res = expr[j] * EXPR.likelihood[[j]] * sum_k * sum_l
      #print(res)
      sample1_G1model_liks <- sample1_G1model_liks + res
    }
    result[i] <- -log(sample1_G1model_liks)
  }
  return(result)
}

# function to generate facData file
fac <- function(exprGene, integrationMeth) {
  pasteFac <- function(xx) paste('[1,',res_expr,']((',paste(xx,sep="",collapse=","),'))',sep="",collapse="")
  tempS_G1 <- matrix(ncol=1+length(IDs_promoter)+length(IDs_body),nrow=length(G1))
  tempS_G1[,1] <- apply(exprGene, 1, pasteFac)
  tempS_G1[,2:(length(IDs_promoter)+1)] <- t(apply(integrationMeth[IDs_promoter, , , drop=F], c(1,2), pasteFac))
  tempS_G1[,(length(IDs_promoter)+2):(1+length(IDs_promoter)+length(IDs_body))] <- t(apply(integrationMeth[IDs_body, , , drop=F], c(1,2), pasteFac))
  return(tempS_G1)
}

for (i in length(workingList)) {
  system(command=paste('mkdir',i,sep=" "))
  system(command=paste('mkdir ./',i,'/G1_model',sep=""))
  system(command=paste('mkdir ./',i,'/G1_model/all',sep=""))
  system(command=paste('mkdir ./',i,'/G2_model',sep=""))
  system(command=paste('mkdir ./',i,'/G2_model/all',sep=""))
  system(command=paste('mkdir ./',i,'/full_model',sep=""))
  system(command=paste('mkdir ./',i,'/null',sep=""))
  system(command=paste('mkdir ./',i,'/null/G1_model',sep=""))
  system(command=paste('mkdir ./',i,'/null/G2_model',sep=""))
  
  gene <- workingList[i]
  IDs_promoter <- eval(parse(text = paste0('prIND$SID$','"',gene,'"')))
  IDs_body <- eval(parse(text = paste0('gbIND$SID$','"',gene,'"')))
  
  #smooth parameters
  smooth_e <- 1/(mean(c(length(G1),length(G2)))/res_expr) # rules of thumb
  smooth_pr <- trunc(1/(mean(c(length(G1),length(G2)))*length(IDs_promoter)/(res_expr*res_expr))) # rules of thumb
  smooth_gb <- trunc(1/(mean(c(length(G1),length(G2)))*length(IDs_body)/(res_expr*res_expr))) # rules of thumb
  sm_param <- c(smooth_e, smooth_gb, smooth_pr)
  #CpG names
  promoter_CpGs <- paste0(template_promoter_CpGs[1:length(IDs_promoter)], '.likelihood')
  geneBody_CpGs <- paste0(template_body_CpGs[1:length(IDs_body)], '.likelihood')
  #start precomputing correct initialization of parameters
  prGeneCpG <- t(apply(integrationMeth[IDs_promoter, , , drop=F],2, function(xx) {xx = apply(xx, 2, geo_mean); xx = xx / sum(xx)}))
  gbGeneCpG <- t(apply(integrationMeth[IDs_body, , , drop=F],2, function(xx) {xx = apply(xx, 2, geo_mean); xx = xx / sum(xx)}))
  exprGene <- integrationExpr[gene,,]
  
  ###########################################################################
  ############################   Tumor model  ###############################
  G1_prior = pot(prGeneCpG[G1,], gbGeneCpG[G1,], exprGene[G1,], sm_param)  
  G1_G1model_mlogliks <- calLike(exprGene[G1, , drop=F], integrationMeth[c(IDs_promoter, IDs_body), G1, , drop=F], G1_prior)

  ###########################################################################
  ###########################  Normal model   ###############################
  G2_prior = pot(prGeneCpG[G2,], gbGeneCpG[G2,], exprGene[G2,], sm_param)  
  G2_G2model_mlogliks <- calLike(exprGene[G2, , drop=F], integrationMeth[c(IDs_promoter, IDs_body), G2, , drop=F], G2_prior)

  ###########################################################################
  ############################   Full model   ###############################
  Full_prior = pot(prGeneCpG[c(G1,G2),], gbGeneCpG[c(G1,G2),], exprGene[c(G1,G2),], sm_param)  
  allData_jointModel_mlogliks <- calLike(exprGene[c(G1,G2), , drop=F], integrationMeth[c(IDs_promoter, IDs_body), c(G1,G2), , drop=F], Full_prior)

  ###########################################################################
  ######################## D calculation ###################################
  D <- 2*(sum(allData_jointModel_mlogliks) - (sum(G1_G1model_mlogliks)+sum(G2_G2model_mlogliks)))
  ###########################################################################
  ################# P val calculation using null distr. #####################
  ###########################################################################
  Ds <- vector(length=nruns,mode="numeric")
  for (run in 1:nruns) {
    cur <- sample(x=1:(length(G2)+length(G1)),size=length(G2),replace=FALSE)
    # G1
    null_G1_prior = pot(prGeneCpG[-cur,], gbGeneCpG[-cur,], exprGene[-cur,], sm_param)  
    G1_G1model_mlogliks <- calLike(exprGene[G1, , drop=F], integrationMeth[c(IDs_promoter, IDs_body), G1, , drop=F], null_G1_prior)
    
    # G2
    null_G2_prior = pot(prGeneCpG[cur,], gbGeneCpG[cur,], exprGene[cur,], sm_param)  
    G2_G2model_mlogliks <- calLike(exprGene[G2, , drop=F], integrationMeth[c(IDs_promoter, IDs_body), G2, , drop=F], null_G2_prior)

    Ds[run] <- 2*(sum(allData_jointModel_mlogliks) - (sum(G1_G1model_mlogliks)+sum(G2_G2model_mlogliks)))
    if(is.na(Ds[run])) break
  }
  #if (sd(Ds, na.rm = T) != 0 & D > 0.1) zscore <- (D - mean(Ds, na.rm = T)) / sd(Ds, na.rm = T) else zscore <- -6
  if (sd(Ds) != 0 & D > 0.1) zscore <- (D - mean(Ds)) / sd(Ds) else zscore <- -6
  pval_zscore <- pnorm(zscore,lower.tail=FALSE)
  ###########################################################################################
  
  eval(parse(text=paste('write.table(x=t(c(pval_zscore,D,mean(Ds),sd(Ds),zscore)), col.names=FALSE, row.names=FALSE, append=TRUE, file="./',i,'.result")',sep="")))
  system(intern=TRUE,command=paste('tar cf ',i,'.tar ',i,sep=""))
}
