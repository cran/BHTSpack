mu.k.u = function(k, ik, z, sigma, mu0){
  nk = sum(ik==k)
  z = z[ik==k]

  if (nk > 0){
    z.mean = mean(z)
    mean.aux = (nk*z.mean+mu0) / (nk+1)
    sigma.aux = sigma[k]/(nk+1)
  }
  else{
    mean.aux = mu0
    sigma.aux = sigma[k]
  }

  return(rnorm(1, mean.aux, sqrt(sigma.aux)))
}
