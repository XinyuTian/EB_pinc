# calculate sample_G1model_mlogliks 
# variables include prior and conditional prop
# prior: expr (1*25) M_exprGB (25*25) M_exprP (25*25) 
# conditional prop: EXPR.likelihood (25), mP.likelihood (6*25), mGB.likelihood (10*25)
# example of calling ---- Rscript calc_lik.R ~/Downloads/samples_FacData1.tab ~/Downloads/factorPotentials.txt 
args <- commandArgs(trailingOnly = TRUE)
file1 <- args[1] # first ID to process
file2 <- args[2] # last consequitve ID to process
### read.Fac
inputFile <- file1
ti <- read.table(inputFile, stringsAsFactors = F, nrows = 1, row.names = 1)
ti_expr <- grep("EXPR", ti)
ti_p <- grep("CpG_P", ti)
ti_gb <- grep("CpG_GB", ti)

tt <- read.table(inputFile, stringsAsFactors = F, skip = 1, row.names = 1)

str2tabs <- function(str, t1=ti_expr, t2=ti_p, t3=ti_gb){
  str = lapply(str, function(xx) sub("\\).*", "", sub(".*\\(", "", xx)) )
  str = lapply(str, function(xx) strsplit(xx, ',') )
  str = matrix(as.numeric(unlist(str)), ncol = 25, byrow = T)
  m1 = str[t1, ]
  m2 = str[t2, ]
  m3 = str[t3, ]
  result = list("EXPR.likelihood"=m1, "mP.likelihood"=m2, "mGB.likelihood"=m3)
  return(result)
}

suppressWarnings(library("plyr"))
facData <- alply(tt, 1, str2tabs) # n.samples * (1+6+10) * 25
### End read.Fac

### read.Potential
inputFile <- file2
temp_files = paste0("M", 1:5, "_temp.txt")
con  <- file(inputFile, open = "r")
flag = rep(F, 5)
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  if(grepl("pot_EXPR.M.GB", oneLine)) flag[1] = T else 
    if(grepl("pot_EXPR.M.P", oneLine)) flag[2] = T else
      if(grepl("pot_EXPR.likelihood", oneLine)) flag[3] = T else
        if(grepl("pot_CpG_P", oneLine)) flag[4] = T else
          if(grepl("pot_CpG_GB", oneLine)) flag[5] = T else
            if(grepl("PC_MAT", oneLine)) flag = rep(F, 5)
  if(any(flag)) eval(parse(text = paste('write.table(', paste('oneLine,file = "',temp_files[which(flag)],'", append = T, quote = F, row.names = F, col.names = F)', sep = ""))))
} 
close(con)

tt = read.table("M1_temp.txt", skip = 2, nrows = 1, stringsAsFactors = F)
M_exprGB = c()
M_exprGB[1] = tt[2]
for(i in 2:25) M_exprGB[i] = read.table("M1_temp.txt", skip = i+1, nrows = 1, stringsAsFactors = F)
M_exprGB = lapply(M_exprGB, function(xx) sub("\\).*", "", sub(".*\\(", "", xx)) )
M_exprGB = lapply(M_exprGB, function(xx) strsplit(xx, ',') )
M_exprGB = matrix(as.numeric(unlist(M_exprGB)), nrow = 25, byrow = T)

tt = read.table("M2_temp.txt", skip = 2, nrows = 1, stringsAsFactors = F)
M_exprP = c()
M_exprP[1] = tt[2]
for(i in 2:25) M_exprP[i] = read.table("M2_temp.txt", skip = i+1, nrows = 1, stringsAsFactors = F)
M_exprP = lapply(M_exprP, function(xx) sub("\\).*", "", sub(".*\\(", "", xx)) )
M_exprP = lapply(M_exprP, function(xx) strsplit(xx, ',') )
M_exprP = matrix(as.numeric(unlist(M_exprP)), nrow = 25, byrow = T)

tt = read.table("M3_temp.txt", skip = 2, nrows = 1, stringsAsFactors = F)
expr = as.numeric(unlist(strsplit(sub("\\).*", "", sub(".*\\(", "", tt[2])), ',')))

mP = list(); i=1
con  <- file("M4_temp.txt", open = "r")
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  if(grepl("POT_MAT:", oneLine)) {mP[i] = oneLine; i = i + 1}
} 
close(con)
mP = lapply(mP, function(xx) sub("\\).*", "", sub(".*\\(", "", xx)) )
mP = lapply(mP, function(xx) strsplit(xx, ',') )
mP = matrix(as.numeric(unlist(mP)), ncol = 25, byrow = T)

mGB = list(); i=1
con  <- file("M5_temp.txt", open = "r")
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  if(grepl("POT_MAT:", oneLine)) {mGB[i] = oneLine; i = i + 1}
} 
close(con)
mGB = lapply(mGB, function(xx) sub("\\).*", "", sub(".*\\(", "", xx)) )
mGB = lapply(mGB, function(xx) strsplit(xx, ',') )
mGB = matrix(as.numeric(unlist(mGB)), ncol = 25, byrow = T)

invisible(file.remove(temp_files))
### End read.Potential


for(i in 1:length(facData)) {
  EXPR.likelihood = facData[[i]][[1]] 
  mP.likelihood = facData[[i]][[2]] 
  mGB.likelihood = facData[[i]][[3]] 
  sample1_G1model_liks = 0
  for (j in seq(25)){
    sum_k = 0 
    for (k in seq(25)) {
      sum_k = sum_k + M_exprP[j,k] * prod(mP.likelihood[, k])
    }
    sum_l = 0
    for (l in seq(25)) {
      sum_l = sum_l + M_exprGB[j,l] * prod(mGB.likelihood[, l])
    }
    res = expr[j] * EXPR.likelihood[[j]] * sum_k * sum_l
    #print(res)
    sample1_G1model_liks <- sample1_G1model_liks + res
  }
  sample1_G1model_mlogliks = -log(sample1_G1model_liks)
  cat(paste0(sample1_G1model_mlogliks, "\n"))
}
