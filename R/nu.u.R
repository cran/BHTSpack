nu.u = function(ind, tau, H){
  n = length(ind)
  nu = rep(0, H)
  nu <- .C("abfun", as.integer(ind-1), as.integer(n), as.double(tau), as.integer(H), as.double(nu), PACKAGE="BHTSpack")[[5]]

  return(nu)
}
