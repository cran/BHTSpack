ptrace = function(res, var, ndisc, nr, nc, type="trace"){
  res = res[["dat.store"]][[var]]

if (type=="trace"){
  if (var == "pk0")
    var = expression(paste({lambda[k]^{(0)}}))
  else if (var == "pk1")
    var = expression(paste({lambda[k]^{(1)}}))
  else if (var == "mu0")
    var = expression(paste({mu[0][k]}))
  else if (var == "mu1")
    var = expression(paste({mu[1][k]}))
  else if (var == "sigma0")
    var = expression(paste({sigma[0][k]^2}))
  else if (var == "sigma1")
    var = expression(paste({sigma[1][k]^2}))
  else
    stop("Unknown var")
}

if (type=="acf"){
  if (var == "pk0")
    var = expression(paste("ACF (", {lambda[k]^{(0)}}, ")", sep=""))
  else if (var == "pk1")
    var = expression(paste("ACF (", {lambda[k]^{(1)}}, ")", sep=""))
  else if (var == "mu0")
    var = expression(paste("ACF (", {mu[0][k]}, ")", sep=""))
  else if (var == "mu1")
    var = expression(paste("ACF (", {mu[1][k]}, ")", sep=""))
  else if (var == "sigma0")
    var = expression(paste("ACF (", {sigma[0][k]^2}, ")", sep=""))
  else if (var == "sigma1")
    var = expression(paste("ACF (", {sigma[1][k]^2}, ")", sep=""))
  else
    stop("Unknown var")
}

  if (!is.null(ncol(res))){
    ncomp = ncol(res)
    if (ncomp > nr*nc)
      stop("Insufficient Number of Plot Cells. Increase nr or nc")

    res = res[-seq(1, ndisc),]
    res = t(apply(res, 1, function(x){sort(x)}))

    layout(matrix(seq(1, nr*nc), nr, nc, byrow=TRUE))
    if (type=="trace")
      invisible(sapply(1:ncomp, function(x){d=res[,x]; plot(d, type="l", xlab="iteration", ylab=var, ylim=range(res), main=paste("k =", x, sep=" "), cex.lab=1);})) 
    if (type=="acf")
      invisible(sapply(1:ncomp, function(x){d=res[,x]; acf(d,lag.max=length(d),ylab=var,main=paste("k =", x, sep=" "), cex.lab=1);})) 
   }
   else{
     res = res[-seq(1, ndisc)]
     plot(res, xlab="iteration", pch=16, ylab=var)
   }
}
