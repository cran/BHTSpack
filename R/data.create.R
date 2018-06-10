data.create = function(N, nr, nc, M, p, s=NULL, covrow=NULL, covcol=NULL, c=0.0001, mat=FALSE){
  set.seed(s)

  nm = N*M

  mscale = c(1, 1.2, 1.4, 1.6)
  vscale = c(1, 1.1, 1.2, 1.3)

  mean00 = 0.1
  mean0 = rep(mean00, 4) * mscale
  var00 = 0.01
  var0 = rep(var00, 4) * vscale

  mean10 = 0.2
  mean1 = rep(mean10, 4) * mscale
  var10 = 0.002
  var1 = rep(var10, 4) * vscale

  hitpr = rep(1/4, 4)

  B = rbinom(n=N*M, size=1, prob=p)
  I = sapply(1:nm, function(x){which(rmultinom(1,1,hitpr)==1)})

  Z = sapply(1:nm, function(x){if(B[x]==0)
                               return(rlnorm(1, lg.mu.sig(mean0[I[x]],var0[I[x]])[["mu"]], lg.mu.sig(mean0[I[x]],var0[I[x]])[["sig"]]))
                             else
                               return(rlnorm(1, lg.mu.sig(mean1[I[x]],var1[I[x]])[["mu"]], lg.mu.sig(mean1[I[x]],var1[I[x]])[["sig"]]))})

  ncs = c(0,cumsum(rep(N, M)))
  Z = lapply(2:length(ncs), function(x){Z[(ncs[x-1]+1):ncs[x]]})
  Z = lapply(Z, function(x){d=matrix(x, nr, nc, byrow=TRUE); rownames(d)=LETTERS[seq(1,nrow(d))]; colnames(d)=seq(1,ncol(d)); return(d);})

  I = lapply(2:length(ncs), function(x){I[(ncs[x-1]+1):ncs[x]]})
  I = lapply(I, function(x){d=matrix(x, nr, nc, byrow=TRUE); rownames(d)=LETTERS[seq(1,nrow(d))]; colnames(d)=seq(1,ncol(d)); return(d);})
  names(I) = paste("Plate", seq(1, length(Z)), sep="")
  I = lapply(I, function(x){d=as.vector(x); names(d)=sapply(colnames(x), function(y){paste(rownames(x), y, sep="")}); return(d);})

  B = lapply(2:length(ncs), function(x){B[(ncs[x-1]+1):ncs[x]]})
  B = lapply(B, function(x){d=matrix(x, nr, nc, byrow=TRUE); rownames(d)=LETTERS[seq(1,nrow(d))]; colnames(d)=seq(1,ncol(d)); return(d);})
  names(B) = paste("Plate", seq(1, length(Z)), sep="")
  B = lapply(B, function(x){d=as.vector(x); names(d)=sapply(colnames(x), function(y){paste(rownames(x), y, sep="")}); return(d);})

  ### incorporate noise
  if (!is.null(covrow) & !is.null(covcol)){
    set.seed(s)

    covrow_chol = chol(covrow)
    covcol_chol = chol(covcol)

    dnoise = lapply(1:M, function(x){matrix(rnorm(nr*nc), nr, nc)})  
    dnoise = lapply(dnoise, function(x){c*covrow_chol%*%x%*%covcol_chol})
    dnoise = lapply(dnoise, function(x){d=as.vector(x); names(d)=sapply(colnames(x), function(y){paste(rownames(x), y, sep="")}); return(d);})

    #add to compound data
    Z = lapply(1:M, function(x){Z[[x]]+dnoise[[x]]})
  }

  names(Z) = paste("Plate", seq(1, length(Z)), sep="")
  if (!mat)
    Z = lapply(Z, function(x){d=as.vector(x); names(d)=sapply(colnames(x), function(y){paste(rownames(x), y, sep="")}); return(d);})

  return(list(Z=Z, I=I, B=B))
}
