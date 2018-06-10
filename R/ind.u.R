ind.u = function(pr){
  n = nrow(pr)
  K = ncol(pr)

  output = rep(0, n) 
  out<-.C("multinomind", as.double(t(pr)), as.integer(n), as.integer(K), as.integer(output), PACKAGE="BHTSpack")
  return(out[[4]])
}
