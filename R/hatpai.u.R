hatpai.u = function(z, hk1, hk0, ph1, ph0, sigma1, sigma0, mu1, mu0, pai, H, n){
  sigma1 = sigma1[hk1]
  mu1 = mu1[hk1]

  sigma0 = sigma0[hk0]
  mu0 = mu0[hk0]

  out <- .Call("hat_pai", z, ph1, ph0, mu1, mu0, sigma1, sigma0, pai, n, H, PACKAGE="BHTSpack")
  return(out)
}

