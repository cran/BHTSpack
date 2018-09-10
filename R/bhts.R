bhts = function(Z, iters, H, K, mu00=NULL, mu10=NULL, a.alpha, b.alpha, a.tau, b.tau, pnorm=FALSE, s=NULL, store=FALSE){
  if (!is.null(s))
    set.seed(s)

  if (pnorm)
    Z = lapply(Z, function(x){(x-mean(x))/sd(x)})

  a0.alpha = a.alpha
  a1.alpha = a.alpha
  b0.alpha = b.alpha
  b1.alpha = b.alpha

  a0.tau = a.tau
  a1.tau = a.tau
  b0.tau = b.tau
  b1.tau = b.tau

  burn = floor(iters/2)
  sampl = 1
  M = length(Z)
  n = unlist(lapply(Z, length))
  z = unlist(Z)

  a = abfun(var(z), 10^-4)[["a"]]
  b = abfun(var(z), 10^-4)[["b"]]

  a0 = a
  a1 = a
  b0 = b
  b1 = b

  if (is.null(mu00) | is.null(mu10)){
    mu00 = mean(z) / 2
    mu10 = 3 * mu00
  }

  # parameters to be saved
  b.s = rep(0, sum(n))
  hatpai.s = rep(0, sum(n))

  # initialize parameters
  a.pai = 0.5*sum(n)
  b.pai = 0.5*sum(n)

  sigma0 = 1/rgamma(n=K, shape=a0, rate=b0)
  mu0 = rnorm(K, mu00, 1)

  sigma1 = 1/rgamma(n=K, shape=a1, rate=b1)
  mu1 = rnorm(K, mu10, 1)

  pai = rbeta(1, a.pai, b.pai)
  hatpai = runif(sum(n))
  b = rbinom(n=length(hatpai), size=1, prob=hatpai)

  ncs = c(0,cumsum(n))
  b.aux = lapply(2:length(ncs), function(x){b[(ncs[x-1]+1):ncs[x]]})

  ih0 = lapply(n, function(x){sample(H, x, replace=TRUE)})
  ih1 = lapply(n, function(x){sample(H, x, replace=TRUE)})

  hk0 = matrix(sample(K, M*H, replace=TRUE), M, H)
  hk1 = matrix(sample(K, M*H, replace=TRUE), M, H)

  alpha0 = rgamma(1, 1, 1)
  alpha1 = rgamma(1, 1, 1)

  tau0 = rgamma(1, 1, 1)
  tau1 = rgamma(1, 1, 1)

  nu.h0 = lapply(ih0, nu.u, alpha0, H)
  nu.h1 = lapply(ih1, nu.u, alpha1, H)

  nu.k0 =  nu.u(hk0, tau0, K)
  nu.k1 =  nu.u(hk1, tau1, K)

  ph0 = lapply(nu.h0, lambda.u)
  ph1 = lapply(nu.h1, lambda.u)

  pk0 = lambda.u(nu.k0)
  pk1 = lambda.u(nu.k1)

  ik0 = unlist(lapply(1:M, function(x){hk0[x,ih0[[x]]]}))
  ik1 = unlist(lapply(1:M, function(x){hk1[x,ih1[[x]]]}))

  h.pr0 = h.pr.u(z[b==0], unlist(ih0)[b==0], mu0, sigma0, pk0, K, H, unlist(lapply(b.aux, function(x){sum(1-x)}))) 
  h.pr1 = h.pr.u(z[b==1], unlist(ih1)[b==1], mu1, sigma1, pk1, K, H, unlist(lapply(b.aux, sum)))

  z.pr0 = z.pr.u(z, as.vector(t(hk0)), mu0, sigma0, unlist(ph0), H, n) 
  z.pr1 = z.pr.u(z, as.vector(t(hk1)), mu1, sigma1, unlist(ph1), H, n) 

  if (store){
    # variables to be stored
    sigma0.store = matrix(rep(0, K*iters), iters, K)
    mu0.store = sigma0.store
    sigma1.store = sigma0.store
    mu1.store = sigma0.store
    pk0.store = sigma0.store
    pk1.store = sigma0.store

    tau1.store = rep(0, iters)
    tau0.store = tau1.store

    nu.k1.store = matrix(rep(0, K*iters), iters, K)
    nu.k0.store = matrix(rep(0, K*iters), iters, K)
  }

  iter = 1
  stop = FALSE
  while (!stop){
    if (iter%%10==0){
      cat("iter=", iter, ", ", sep="")

      if (iter%%100==0)
        cat("\n")
    }
 
    if (store){
      sigma0.store[iter,] = sigma0
      sigma1.store[iter,] = sigma1

      mu0.store[iter,] = mu0
      mu1.store[iter,] = mu1

      pk0.store[iter,] = pk0
      pk1.store[iter,] = pk1
      tau1.store[iter] = tau1
      tau0.store[iter] = tau0

      nu.k1.store[iter,] = nu.k1
      nu.k0.store[iter,] = nu.k0
    }

    # update pai
    pai = pai.u(b, a.pai, b.pai)

    # update hatpai
    hatpai = hatpai.u(z, as.vector(t(hk1)), as.vector(t(hk0)), unlist(ph1), unlist(ph0), sigma1, sigma0, mu1, mu0, pai, H, n)

    # update b
    b = b.u(hatpai)

    # update sigma
    sigma0 = sapply(1:K, sig.k.u, ik0[b==0], z[b==0], mu00, a0, b0)
    sigma1 = sapply(1:K, sig.k.u, ik1[b==1], z[b==1], mu10, a1, b1)

    # update mu
    mu0 = sapply(1:K, mu.k.u, ik0[b==0], z[b==0], sigma0, mu00)
    mu1 = sapply(1:K, mu.k.u, ik1[b==1], z[b==1], sigma1, mu10)

    # update z.pr
    z.pr0 = z.pr.u(z, as.vector(t(hk0)), mu0, sigma0, unlist(ph0), H, n)
    z.pr1 = z.pr.u(z, as.vector(t(hk1)), mu1, sigma1, unlist(ph1), H, n)

    # update ih
    ih0 = ind.u(z.pr0)
    ih1 = ind.u(z.pr1)

    ncs = c(0,cumsum(n))
    b.aux = lapply(2:length(ncs), function(x){b[(ncs[x-1]+1):ncs[x]]})

    # update h.pr
    if (sum(b==0)!=0)
      h.pr0 = h.pr.u(z[b==0], ih0[b==0], mu0, sigma0, pk0, K, H, unlist(lapply(b.aux, function(x){sum(1-x)}))) 
    else
      h.pr0 = matrix(rep(1/K, M*H*K), M*H, K)

    if (sum(b==1)!=0)
      h.pr1 = h.pr.u(z[b==1], ih1[b==1], mu1, sigma1, pk1, K, H, unlist(lapply(b.aux, sum)))
    else
      h.pr1 = matrix(rep(1/K, M*H*K), M*H, K)
    
    # reorganize ih into list
    ih0 = lapply(2:length(ncs), function(x){ih0[(ncs[x-1]+1):ncs[x]]})
    ih1 = lapply(2:length(ncs), function(x){ih1[(ncs[x-1]+1):ncs[x]]})

    # update hk
    hk0 = matrix(ind.u(h.pr0), M, H, byrow=TRUE)
    hk1 = matrix(ind.u(h.pr1), M, H, byrow=TRUE)

    # update ik
    ik0 = unlist(lapply(1:M, function(x){hk0[x,ih0[[x]]]}))
    ik1 = unlist(lapply(1:M, function(x){hk1[x,ih1[[x]]]}))

    # update nu.h
    nu.h0 = lapply(1:M, function(x){nu.u(ih0[[x]][b.aux[[x]]==0], alpha0, H)})
    nu.h1 = lapply(1:M, function(x){nu.u(ih1[[x]][b.aux[[x]]==1], alpha1, H)})

    # update nu.k
    nu.k0 = nu.u(hk0, tau0, K)
    nu.k1 = nu.u(hk1, tau1, K)

    # update alpha
    alpha0 = alpha.u(nu.h0, a0.alpha, b0.alpha, H)
    alpha1 = alpha.u(nu.h1, a1.alpha, b1.alpha, H)

    # update tau
    tau0 = tau.u(nu.k0, a0.tau, b0.tau)
    tau1 = tau.u(nu.k1, a1.tau, b1.tau)

    # update ph
    ph0 = lapply(nu.h0, lambda.u)
    ph1 = lapply(nu.h1, lambda.u)

    # update pk
    pk0 = lambda.u(nu.k0)
    pk1 = lambda.u(nu.k1)  

    if (iter > burn & iter%%sampl==0){
      b.s = b.s + b
      hatpai.s = hatpai.s + hatpai
    }

    if (iter == iters)
      stop = TRUE
    else
      iter = iter + 1
  }

  sampl.num = as.integer((iters-burn)/sampl)

  hatpai = hatpai.s / sampl.num

  ncs = c(0,cumsum(n))

  hatpai = lapply(2:length(ncs), function(x){hatpai[(ncs[x-1]+1):ncs[x]]})
  hatpai = lapply(1:length(n), function(x){v=hatpai[[x]]; names(v)=names(Z[[x]]); return(v);})
  names(hatpai) = names(Z)

  if (store){
    dat.store = list(sigma0=sigma0.store, 
                     sigma1 = sigma1.store, 
                     mu0 = mu0.store,
                     mu1 = mu1.store,
                     pk0 = pk0.store,
                     pk1 = pk1.store,
                     tau1 = tau1.store,
                     tau0 = tau0.store,
                     nu.k1 = nu.k1.store,
                     nu.k0 = nu.k0.store)
    return(list(hatpai=hatpai, dat.store=dat.store))
  }
  else
    return(list(hatpai=hatpai))
}
