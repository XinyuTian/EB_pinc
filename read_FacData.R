setwd("~/Downloads/")
inputFile <- "./samples_FacData1.tab"
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

library(plyr)
facData <- alply(tt, 1, str2tabs) # n.samples * (1+6+10) * 25
