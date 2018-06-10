alpha.u = function(nu, a0, b0, H){
  M = length(nu)
  a = a0 + M*(H-1)
  b = b0 - sum(unlist(lapply(nu, function(x){log(1-x[-H])})))

  return(rgamma(1, a, b))
}
