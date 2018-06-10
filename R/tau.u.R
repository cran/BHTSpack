tau.u = function(nu, a0, b0){
  K = length(nu)
  a = a0 + K - 1
  b = b0 - sum(log(1-nu[-K]))

  return(rgamma(1, a, b))
}
