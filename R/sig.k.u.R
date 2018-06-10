sig.k.u = function(k, ik, z, mu0, a0, b0){
  nk = sum(ik==k)
  z = z[ik==k]

  if (nk > 0){
    z.mean = mean(z)
    a = a0+nk
    b = b0 + 0.5*nk/(nk+1)*((mu0-z.mean)^2) + 0.5*sum((z-z.mean)^2)
  }
  else{
    a = a0
    b = b0
  }

  return(1/rgamma(n=1, shape=a, rate=b))
}
