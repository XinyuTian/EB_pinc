## use edgeR to get shrinked variance for sequencing data
library(edgeR)
seq.data <- read.csv('gene_counts.csv', row.names = 1) # 23373
condition <- c(rep(1,3), rep(2, 3), rep(0,3))
design <- model.matrix(~condition)
dge <- DGEList(counts=seq.data, group = condition)
# filter
seq.data.use <- seq.data[rowSums(seq.data)>10, ] # 16638
# CPM
logcpm <- cpm(seq.data.use, prior.count=0.1, log=TRUE)

############  -----------   0 vs 5  -----------   ############
condition1 <- c(rep(1,3), rep(0,3))
design1 <- model.matrix(~condition1)
## test - two groups
dge1 <- DGEList(counts=seq.data.use[c(1:3,7:9)], group=condition1)
dge1 <- calcNormFactors(dge1) # normalization
dge1 <- estimateCommonDisp(dge1)
dge1 <- estimateTagwiseDisp(dge1)
et1 <- exactTest(dge1) 
topTags(et1)
## GLM - two groups
## dge1$trended.dispersion is the shrinked dispersion
dge1 <- DGEList(counts=seq.data.use[c(1:3,7:9)], group=condition1)
dge1 <- estimateGLMCommonDisp(dge1, design1)
dge1 <- estimateGLMTrendedDisp(dge1, design1)
dge1 <- estimateGLMTagwiseDisp(dge1, design1)
et1 <- exactTest(dge1) 
topTags(et1)

############  -----------   0 vs 5 vs 10 -----------   ############
## GLM - three groups
condition2 <- c(rep(1,3), rep(2,3), rep(0,3))
design2 <- model.matrix(~condition)
dge2 <- DGEList(counts=seq.data.use, group=condition2)
dge2 <- estimateGLMCommonDisp(dge2, design2)
dge2 <- estimateGLMTrendedDisp(dge2, design2)
dge2 <- estimateGLMTagwiseDisp(dge2, design2)
et2 <- exactTest(dge2) 
topTags(et2)

###############     standard procedure    #############
cutoffFDR <- 0.05
p.adj.all <- p.adjust(et$table$PValue, method="BH")
ind.FDR <- p.adj.all < cutoffFDR   ## FDR < 0.1
lfc = et$table$logFC[ind.FDR]
ind1 = lfc > 1
ind2 = lfc < -1 
ind.FC = ind1 | ind2 ## fold change > 2 or < 0.5 
geneNamesSelecte <- geneNames[ind.FDR][ind.FC]
