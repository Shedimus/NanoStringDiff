rnegbinom <- function(n, mu = 1, phi = 0.01) {
    rpois(n, rgamma(n, shape = 1/phi, scale = mu * phi))
}


est.dispersion <- function(Y, Y_nph, lamda_i, c, d) {
    
    nsamples = ncol(Y)
    
    mu = rowMeans(Y_nph)
    muY = sweep(matrix(rep(mu, nsamples), ncol = nsamples), 2, c * d, FUN = "*")
    
    
    get.phihat <- function(dat) {
        y = dat[1:nsamples]
        Ey = dat[nsamples + (1:nsamples)]
        
        obj = function(phi) {
            alpha = 1/phi
            tmp1 = 1/(1 + Ey * phi)
            tmp2 = 1 - tmp1
            tmp2[tmp2==0] = 1e-08
            tmp3.p = fun5(cbind(matrix(y, ncol = 1), matrix(lamda_i, ncol = 1), 
                                matrix(tmp2, ncol = 1),rep(alpha,length(y))))
            -tmp3.p
        }
        return(optimize(obj, interval = c(1e-08, 1e+08))$minimum)
    }
    
    
    phi.g = apply(cbind(matrix(Y, ncol = nsamples), 
                        matrix(muY, ncol = nsamples)), 1, get.phihat)
    
    lphi = log(phi.g)
    
    list(phi = phi.g, eta = lphi)
}


compute.baseSigma <- function(phi0, Y, muY, nsamples) {
    set.seed(123)
    nsim = 5
    res = rep(0, nsim)
    c = rep(1, nsamples)
    d = rep(1, nsamples)
    lamda = rep(20, nsamples)
    for (i in 1:nsim) {
        Ysim <- matrix(rnegbinom(length(Y), muY, phi = phi0), ncol = nsamples)
        ## assume 20 is background noise
        lag.sim = est.dispersion(Ysim + 20, Ysim, lamda, c, d)$eta
        sigma = IQR(lag.sim, na.rm = TRUE)/1.349
        res[i] = sigma^2
    }
    mean(res)
} 
