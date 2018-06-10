h.pr.u = function(z, ih, mu, sigma, pk, K, H, n){
    M = length(n)

    out <- .Call("stick_multnorm_h", z, ih-1, pk, sigma, mu, n, H, PACKAGE="BHTSpack")
    return(matrix(out, H*M, K, byrow=TRUE))
}
