fdr.r = function(r, hatpai, fdr){
  res = sapply(hatpai, function(x){ifelse(x>r,1,0)})
  res = sum(res*(1-hatpai)) / sum(res)

  ifelse(!is.nan(res), return(res-fdr), return(1))
}
