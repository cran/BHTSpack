lg.mu.sig = function(m, v){

  mu = log(m^2/sqrt(m^2+v))
  sig = sqrt(log((m^2+v)/m^2))

  return(list(mu=mu, sig=sig))
}
