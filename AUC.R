calculate.AUC.pROC = function(y, prob){
  rprob = rank(prob)
  n1 = sum(y)
  n0 = length(y) - n1
  u = sum(rprob[y == 1]) - n1 * (n1 + 1)/2
  auc=u/(n1 * n0)
  return(auc)
}
