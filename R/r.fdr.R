r.fdr = function(res, fdr=0.05){
  hatpai = unlist(res[["hatpai"]])

  if (fdr < 1-max(hatpai))
    stop("minimum achievable FDR: ", 1-max(hatpai))

  res = data.frame(ID=names(hatpai), hatpai)
  ind = sort(res[["hatpai"]], decreasing=TRUE, index.return=TRUE)[["ix"]]
  res = res[ind,]
  rownames(res) = NULL

  r = uniroot(fdr.r, hatpai=res[["hatpai"]], fdr=fdr, interval=c(0,max(hatpai)-10^-4))[["root"]]
  res = subset(res, hatpai>r)

  return(list(res=res, r=r))
}
