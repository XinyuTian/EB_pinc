## calculate prediction accuracy
## input 1.predicted n1=27, n2=243

tt = read.table('1.predicted', header=T)
tt = tt[,2] > 0.5
t1 = read.table('~/Downloads/pincage/1.predicted', header=T)
t1 = t1[,2] > 0.5
t0 = c(rep(F,27), rep(T,243))
table(t0,tt)
table(t0,t1)