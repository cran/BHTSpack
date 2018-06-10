z.pr.u = function(z, hk, mu, sigma, ph, H, n){
    sigma = sigma[hk]
    mu = mu[hk]

    out <- .Call("stick_multnorm_z", z, ph, sigma, mu, n, H, PACKAGE="BHTSpack")
    return(matrix(out, sum(n), H, byrow=TRUE))
}
