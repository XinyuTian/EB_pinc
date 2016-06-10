setwd("~/Downloads/")
inputFile <- "./factorPotentials.txt"
temp_files = paste0("M", 1:5, "_temp.txt")
file.remove(temp_files)
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

file.remove(temp_files)
