lambda.u = function(nu){
  H = length(nu)
  ph = rep(0, H)
  ph <- .C("lambda", as.double(nu), as.integer(H), as.double(ph), PACKAGE="BHTSpack")[[3]]

  return(ph)
}
