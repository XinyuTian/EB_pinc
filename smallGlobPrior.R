## GOES WITH CALC_LIK.R 

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
  
  prior_pr <- apply(prGeneCpG,2,mean)
  prior_gb <- apply(gbGeneCpG,2,mean)
  prior_expr <- kernsm(apply(exprGene,2,mean),h=smooth_e)
  prior_expr <- prior_expr@yhat/sum(prior_expr@yhat)
  
  string <- paste(prior_pr,collapse=",")
  promoterPots <- paste("\nNAME:\t\tpot_",promoter_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n",sep="",collapse="")
  string <- paste(prior_gb,collapse=",")
  geneBodyPots <- paste("\nNAME:\t\tpot_",geneBody_CpGs,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n",sep="",collapse="")
  string <- paste(prior_expr,collapse=",")
  expr.pots <- paste("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\nNAME:\t\tpot_EXPR.prior\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,",res_expr,"]((",string,"))\nPC_MAT:\t\t[1,",res_expr,"]((",paste(rep(1,res_expr),collapse=","),"))\n\n",sep="",collapse="")
  
  result <- tensor_product(gbGeneCpG,exprGene,smooth_h=smooth_gb)
  expr.m <- paste("NAME:\t\tpot_EXPR.M.GB\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_expr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_expr),collapse=","),"](())\n\n",sep="",collapse="")
  
  result <- tensor_product(prGeneCpG,exprGene,smooth_h=smooth_pr)
  expr.m <- c(expr.m,paste("NAME:\t\tpot_EXPR.M.P\nTYPE:\t\trowNorm\nPOT_MAT:\t\t[",paste(c(res_expr,res_expr),collapse=","),"]((",paste(apply(result,1,paste,collapse=","),collapse="),\n\t\t\t("),"))\nPC_MAT:\t\t[",paste(c(res_expr,res_expr),collapse=","),"](())\n\n",sep="",collapse=""))
  
  return(list = c(expr.m,expr.pots,promoterPots,geneBodyPots))
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
  potentials <- file(paste("./",i,"/G1_model/all/factorPotentials.txt",sep=""),"w")
  cat(G1_prior,file=potentials)
  close(potentials)
  
  tempS_G1 <- fac(exprGene[G1, , drop=F], integrationMeth[c(IDs_promoter, IDs_body), G1, , drop=F])
  rownames(tempS_G1) <- G1
  colnames(tempS_G1) <- c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs)
  eval(parse(text = paste('write.table(', paste('tempS_G1,file ="./',i,'/G1_model/all/G1_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
  
  string<-system(intern=TRUE,command=paste('Rscript ~/Dropbox/My\\ ','R\\ ','Code/EB_pinc/calc_lik.R ./',i,'/G1_model/all/G1_FacData.tab ./',i,'/G1_model/all/factorPotentials.txt',sep=""))
  G1_G1model_mlogliks <- as.numeric(string)
  ###########################################################################
  ###########################  Normal model   ###############################
  G2_prior = pot(prGeneCpG[G2,], gbGeneCpG[G2,], exprGene[G2,], sm_param)  
  potentials <- file(paste("./",i,"/G2_model/all/factorPotentials.txt",sep=""),"w")
  cat(G2_prior,file=potentials)
  close(potentials)
  
  tempS_G2 <- fac(exprGene[G2, , drop=F], integrationMeth[c(IDs_promoter, IDs_body), G2, , drop=F])
  rownames(tempS_G2) <- G2
  colnames(tempS_G2) <- c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs)
  eval(parse(text = paste('write.table(', paste('tempS_G2,file ="./',i,'/G2_model/all/G2_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
  
  string<-system(intern=TRUE,command=paste('Rscript ~/Dropbox/My\\ ','R\\ ','Code/EB_pinc/calc_lik.R ./',i,'/G2_model/all/G2_FacData.tab ./',i,'/G2_model/all/factorPotentials.txt',sep=""))
  G2_G2model_mlogliks <- as.numeric(string)
  ###########################################################################
  ############################   Full model   ###############################
  Full_prior = pot(prGeneCpG[c(G1,G2),], gbGeneCpG[c(G1,G2),], exprGene[c(G1,G2),], sm_param)  
  potentials <- file(paste("./",i,"/full_model/factorPotentials.txt",sep=""),"w")
  cat(Full_prior,file=potentials)
  close(potentials)
  
  tempFac <- rbind(tempS_G1,tempS_G2)
  rownames(tempFac) <- c(G1,G2)
  colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs,geneBody_CpGs)
  eval(parse(text = paste('write.table(', paste('tempFac,file = "./',i,'/full_model/full_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
  
  # query the full model with T and AN samples
  string<-system(intern=TRUE,command=paste('Rscript ~/Dropbox/My\\ ','R\\ ','Code/EB_pinc/calc_lik.R ./',i,'/full_model/full_FacData.tab ./',i,'/full_model/factorPotentials.txt',sep=""))
  allData_jointModel_mlogliks <- as.numeric(string)
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
    potentials <- file(paste("./",i,"/null/G1_model/factorPotentials.txt",sep=""),"w")
    cat(null_G1_prior,file=potentials)
    close(potentials)
    
    # query
    string<-system(intern=TRUE,command=paste('Rscript ~/Dropbox/My\\ ','R\\ ','Code/EB_pinc/calc_lik.R ./',i,'/G1_model/all/G1_FacData.tab ./',i,'/null/G1_model/factorPotentials.txt',sep=""))
    G1_G1model_mlogliks <- as.numeric(string)

    # G2
    null_G2_prior = pot(prGeneCpG[cur,], gbGeneCpG[cur,], exprGene[cur,], sm_param)  
    potentials <- file(paste("./",i,"/null/G2_model/factorPotentials.txt",sep=""),"w")
    cat(null_G2_prior,file=potentials)
    close(potentials)
    
    # query
    string<-system(intern=TRUE,command=paste('Rscript ~/Dropbox/My\\ ','R\\ ','Code/EB_pinc/calc_lik.R ./',i,'/G2_model/all/G2_FacData.tab ./',i,'/null/G2_model/factorPotentials.txt',sep=""))
    G2_G2model_mlogliks <- as.numeric(string)

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
