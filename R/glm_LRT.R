glm.LRT <- function(NanoStringData, design.full, Beta = ncol(design.full), contrast = NULL) {
    
    c = positiveFactor(NanoStringData)
    d = housekeepingFactor(NanoStringData)
    k = c * d
    lamda_i = negativeFactor(NanoStringData)
    
    
    if (length(k) == 0) {
        stop("Before calling function glm.LRT, should get normalization factors \n             
             first using function estNormalizationFactors")
    }
    
    
    Y = exprs(NanoStringData)
    Y_n = sweep(Y, 2, lamda_i, FUN = "-")
    Y_nph = sweep(Y_n, 2, k, FUN = "/")
    Y_nph[Y_nph <= 0] = 0.1
    
    X.full = design.full
    nsamples = ncol(Y)
    Beta.names = colnames(design.full)
    
    
    result.full = glmfit.full(NanoStringData, design.full)
    Beta.full = result.full$Beta.full
    U.full = result.full$mean.full
    phi.hat = result.full$dispersion
    df.full = result.full$df.full
    m0 = result.full$m0
    sigma = result.full$sigma
    
    V.full = sweep(U.full, 2, k, FUN = "*")
    
    
    ## Make reduced design matrix.
    ## Here we borrow the idea from paskage edgeR to make reduced design matrix
    
    
    
    
    if (is.null(contrast)) {
        if (length(Beta) > 1) 
            Beta = unique(Beta)
        
        if (is.character(Beta)) {
            check.Beta = Beta %in% Beta.names
            if (any(!check.Beta)) 
                stop("The name(s) of Beta arguments do not match the \n
                     name(s) of the design matrix.")
            Beta = match(Beta, Beta.names)
        }
        
        logFC = Beta.full[, Beta, drop = FALSE]/log(2)
    } else {
        contrast = as.matrix(contrast)
        qrc = qr(contrast)
        ncontrasts = qrc$rank
        
        if (ncontrasts == 0) 
            stop("Need at least one none zero contrast")
        
        Beta = 1:ncontrasts
        
        if (ncontrasts < ncol(contrast)) 
            contrast = contrast[, qrc$pivot[Beta]]
        
        logFC = drop((Beta.full %*% contrast)/log(2))
        
        
        Dvec = rep.int(1, nsamples)
        Dvec[Beta] = diag(qrc$qr)[Beta]
        Q = qr.Q(qrc, complete = TRUE, Dvec = Dvec)
        design.full = design.full %*% Q
    }
    
    
    design.reduce = design.full[, -Beta, drop = FALSE]
    
    
    if (ncol(design.reduce) == 1) {
        
        result.reduce = glmfit.OneGroup(NanoStringData, m0, sigma, phi.hat)
        
    } else {
        
        result.reduce = glmfit.reduce(NanoStringData, design.reduce, m0, sigma, phi.hat)
        
    }
    
    
    Beta.reduce = result.reduce$Beta.reduce
    U.reduce = result.reduce$mean.reduce
    df.reduce = result.reduce$df.reduce
    
    V.reduce = sweep(U.reduce, 2, k, FUN = "*")
    
    
    get.loglikelihood <- function(dat) {
        
        y = dat[1:nsamples]
        Ey = dat[nsamples + (1:nsamples)]
        phi = dat[2 * nsamples + 1]
        
        
        
        alpha = 1/phi
        tmp1 = 1/(1 + Ey * phi)
        tmp2 = 1 - tmp1
        tmp2[tmp2==0] = 1e-08
        
        
        item1 = function(yy) {
            y_gi = yy[1]
            lamda_gi = yy[2]
            tmp2_gi = yy[3]
            t = c(0:y_gi)
            tmp33.t = exp(lgamma(t + alpha) + (y_gi - t) * log(lamda_gi) + t * log(tmp2_gi) - 
                              lfactorial(t) - lfactorial(y_gi - t))
            tmp33.tt = log(max(sum(tmp33.t), 1e-08))
        }
        
        
        tmp3 = apply(cbind(matrix(y, ncol = 1), matrix(lamda_i, ncol = 1), matrix(tmp2, 
                                                                                  ncol = 1)), 1, item1)
        
        sum(tmp3) - nsamples * lgamma(alpha) + alpha * sum(log(tmp1)) - sum(lamda_i)
    }
    
    
    ## compute likelihood under null
    tmpl1 = cbind(Y, V.reduce, phi.hat)
    l0 = apply(tmpl1, 1, get.loglikelihood)
    
    
    ## compute likelihood under alternative
    
    tmpl2 = cbind(Y, V.full, phi.hat)
    la = apply(tmpl2, 1, get.loglikelihood)
    
    
    lr = -2 * (l0 - la)
    
    lr[which(lr <= 0)] = 0
    
    df = df.full - df.reduce
    pval = 1 - pchisq(lr, df = df)
    
    qval = p.adjust(pval, method = "BH")
    
    if (length(Beta) == 1) 
        logFC <- as.vector(logFC)
    
    table = data.frame(logFC = logFC, lr = lr, pvalue = pval, qvalue = qval)
    
    
    list(table = table, dispersion = phi.hat, log.dispersion = log(phi.hat), design.full = X.full, 
         design.reduce = design.reduce, Beta.full = Beta.full, mean.full = U.full, 
         Beta.reduce = Beta.reduce, mean.reduce = U.reduce, m0 = m0, sigma = sigma)
    
    }


glmfit.full <- function(NanoStringData, design.full) {
    
    
    c = positiveFactor(NanoStringData)
    d = housekeepingFactor(NanoStringData)
    k = c * d
    lamda_i = negativeFactor(NanoStringData)
    Y = exprs(NanoStringData)
    Y_n = sweep(Y, 2, lamda_i, FUN = "-")
    Y_nph = sweep(Y_n, 2, k, FUN = "/")
    Y_nph[Y_nph <= 0] = 0.1
    
    
    
    
    nsamples = ncol(Y)
    ngenes = nrow(Y)
    nbeta = ncol(design.full)  # number of full parameters (Beta)
    
    # Beta matrix from linear model, starting value for Betas in optim
    Blm = matrix(NA, ngenes, nbeta)
    
    for (i in 1:ngenes) {
        
        model = lm(log(Y_nph[i, ]) ~ 0 + design.full)
        Blm[i, ] = model$coefficients
        
    }
    
    U = exp(Blm %*% t(design.full))
    V = sweep(U, 2, k, FUN = "*")
    
    phi.g = est.dispersion(Y, Y_nph, lamda_i, c, d)$phi
    
    
    ii = rowMins(Y) > max(negativeControl(NanoStringData))
    
    
    l = length(which(ii == TRUE))
    if (l > 0) {
        
        phi.g0 = phi.g[ii]
        lphi.g0 = log(phi.g0)
        m0 = median(lphi.g0, na.rm = TRUE)
        sigma2.mar = (IQR(lphi.g0, na.rm = TRUE)/1.349)^2
        # Here we borrow the idea to compute the base sigma for DSS The function
        # compute.baseSigma borrow the idea from Hao Wu's function
        # compute.baseSigma.nontrend in DSS Package
        sigma2.base = compute.baseSigma(exp(m0), Y[ii, ], V[ii, ], nsamples)
        sigma = sqrt(max(sigma2.mar - sigma2.base, 0.01))
    } else {
        cat("There is no data satisied that min of endo great than max 
            of negative control ", "\n")
        m0 = -2
        sigma = 1
        lphi.g0 = 10
    }
    
    
    max.phi = max(lphi.g0, 10, na.rm = TRUE)
    max.mean = max(rowMeans(Y_nph))
    
    
    
    
    get.phi <- function(dat) {
        y = dat[1:nsamples]
        Ey = dat[nsamples + (1:nsamples)]
        
        obj = function(phi) {
            alpha = 1/phi
            tmp1 = 1/(1 + Ey * phi)
            tmp2 = 1 - tmp1
            tmp2[tmp2==0] = 1e-08
            
            item1 = function(yy) {
                y_gi = yy[1]
                lamda_gi = yy[2]
                tmp2_gi = yy[3]
                t = c(0:y_gi)
                com = matrix(700, length(t), 1)
                
                
                tmp33.t = exp(rowMins(cbind(lgamma(t + alpha) + (y_gi - t) * log(lamda_gi) + 
                                                t * log(tmp2_gi) - lfactorial(t) - lfactorial(y_gi - t), com)))
                tmp33.tt = log(max(sum(tmp33.t), 1e-08))
                
            }
            tmp3 = apply(cbind(matrix(y, ncol = 1), matrix(lamda_i, ncol = 1), matrix(tmp2, 
                                                                                      ncol = 1)), 1, item1)
            
            
            -(sum(tmp3) - nsamples * lgamma(alpha) + alpha * sum(log(tmp1)) - ((log(phi) - 
                                                                                    m0)^2)/(2 * (sigma^2)) - log(sigma) - sum(lamda_i))
        }
        return(optimize(obj, interval = c(0.005, max.phi))$minimum)
    }
    
    get.beta.full <- function(dat) {
        
        n = nsamples
        y = dat[1:n]
        phi = dat[n + 1]
        Bstart = dat[(n + 2):(n + 1 + nbeta)]
        
        obj = function(beta) {
            
            
            alpha = 1/phi
            xb = beta %*% t(design.full)
            xb[xb > 700] = 700       ## control upper band for exp operation
            tmp1 = 1/(1 + exp(xb) * k * phi)
            tmp2 = 1 - tmp1
            tmp2[tmp2==0] = 1e-08
            
            
            item1 = function(yy) {
                
                
                y_gi = yy[1]
                lamda_gi = yy[2]
                tmp2_gi = yy[3]
                t = c(0:y_gi)
                com = matrix(700, length(t), 1)
                
                
                tmp33.t = exp(rowMins(cbind(lgamma(t + alpha) + (y_gi - t) * log(lamda_gi) + 
                                                t * log(tmp2_gi) - lfactorial(t) - lfactorial(y_gi - t), com)))
                tmp33.tt = log(max(sum(tmp33.t), 1e-08))
            }
            
            tmp3 = apply(cbind(matrix(y, ncol = 1), matrix(lamda_i, ncol = 1), matrix(tmp2, 
                                                                                      ncol = 1)), 1, item1)
            
            
            
            -(sum(tmp3) - n * lgamma(alpha) + alpha * sum(log(tmp1)) - ((log(phi) - 
                                                                             m0)^2)/(2 * (sigma^2)) - log(sigma) - sum(lamda_i))
        }
        
        
        return(optim(Bstart, obj)$par)
    }
    
    
    
    id = c(1:ngenes)
    Beta.full = matrix(0, ngenes, nbeta)
    phi.full = rep(0, ngenes)
    
    phi.s = apply(cbind(matrix(Y, ncol = nsamples), matrix(V, ncol = nsamples)), 
                  1, get.phi)
    B.s = Blm
    Y.t = Y
    
    con11 = 1
    con21 = 1
    
    j = 0
    while ((con11 >= 0.5 | con21 >= 0.001) & j < 50) {
        j = j + 1
        Beta = apply(cbind(matrix(Y.t, ncol = nsamples), matrix(phi.s, ncol = 1), 
                           matrix(B.s, ncol = nbeta)), 1, get.beta.full)
        
        xb = t(design.full %*% Beta)
        xb[xb > 700] = 700         ## control the upper band of exp operation
        U.t = exp(xb)
        V.t = sweep(U.t, 2, k, FUN = "*")
        
        phi.t = apply(cbind(matrix(Y.t, ncol = nsamples), matrix(V.t, ncol = nsamples)), 
                      1, get.phi)
        con1 = rowMaxs(abs((B.s - t(Beta))/t(Beta)))
        con2 = abs((phi.s - phi.t)/phi.s)
        con11 = max(con1)
        con21 = max(con2)
        
        
        idx = which(con1 < 0.5 & con2 < 0.001)
        
        if (!length(idx) == 0) {
            Beta.full[id[idx], ] = t(Beta)[idx, ]
            phi.full[id[idx]] = phi.t[idx]
        }
        
        phi.s = phi.t
        B.s = t(Beta)
        
        if (!length(idx) == 0 & !length(idx) == length(id)) {
            Y.t = Y.t[-idx, ]
            phi.s = phi.s[-idx]
            B.s = B.s[-idx, ]
            id = id[-idx]
        }
        
        # print(j)
        
    }
    
    
    if (j == 50 & !length(idx) == length(id)) {
        
        if(length(idx) == 0) {
            Beta.full[id, ] = t(Beta)
            phi.full[id] = phi.t
        } else{
            Beta.full[id, ] = t(Beta)[-idx,]
            phi.full[id] = phi.t[-idx]
        }
        
        
    }
    
    U.full = exp(Beta.full %*% t(design.full))
    V.full = sweep(U.full, 2, k, FUN = "*")
    
    
    eta = log(phi.full)
    
    
    list(Beta.full = Beta.full, design = design.full, dispersion = phi.full, log.dispersion = eta, 
         m0 = m0, sigma = sigma, df.full = nbeta, mean.full = U.full, nineration = j)
    
    
}



glmfit.reduce <- function(NanoStringData, design.reduce, m0, sigma, phi) {
    
    c = positiveFactor(NanoStringData)
    d = housekeepingFactor(NanoStringData)
    k = c * d
    lamda_i = negativeFactor(NanoStringData)
    Y = exprs(NanoStringData)
    Y_n = sweep(Y, 2, lamda_i, FUN = "-")
    Y_nph = sweep(Y_n, 2, k, FUN = "/")
    Y_nph[Y_nph <= 0] = 0.1
    
    
    
    nsamples = ncol(Y)
    ngenes = nrow(Y)
    nbeta = ncol(design.reduce)  # number of parameters (Beta)
    
    # Beta matrix from linear model, starting value for Betas in optim
    Blm = matrix(NA, ngenes, nbeta)
    
    for (i in 1:ngenes) {
        
        model = lm(log(Y_nph[i, ]) ~ 0 + design.reduce)
        Blm[i, ] = model$coefficients
        
    }
    
    get.beta.reduce <- function(dat) {
        
        n = nsamples
        y = dat[1:n]
        phi = dat[n + 1]
        Bstart = dat[(n + 2):(n + 1 + nbeta)]
        
        obj = function(beta) {
            
            
            alpha = 1/phi
            xb = beta %*% t(design.reduce)
            xb[xb > 700] = 700        ## control upper band for exp operation
            tmp1 = 1/(1 + exp(xb) * k * phi)
            tmp2 = 1 - tmp1
            tmp2[tmp2==0] = 1e-08
            
            
            item1 = function(yy) {
                
                
                y_gi = yy[1]
                lamda_gi = yy[2]
                tmp2_gi = yy[3]
                t = c(0:y_gi)
                com = matrix(700, length(t), 1)
                
                
                tmp33.t = exp(rowMins(cbind(lgamma(t + alpha) + (y_gi - t) * log(lamda_gi) + 
                                                t * log(tmp2_gi) - lfactorial(t) - lfactorial(y_gi - t), com)))
                tmp33.tt = log(max(sum(tmp33.t), 1e-08))
            }
            
            tmp3 = apply(cbind(matrix(y, ncol = 1), matrix(lamda_i, ncol = 1), matrix(tmp2, 
                                                                                      ncol = 1)), 1, item1)
            
            
            -(sum(tmp3) - n * lgamma(alpha) + alpha * sum(log(tmp1)) - ((log(phi) - 
                                                                             m0)^2)/(2 * (sigma^2)) - log(sigma) - sum(lamda_i))
        }
        
        
        return(optim(Bstart, obj)$par)
    }
    
    
    Beta.reduce = apply(cbind(matrix(Y, ncol = nsamples), matrix(phi, ncol = 1), 
                              matrix(Blm, ncol = nbeta)), 1, get.beta.reduce)
    U.reduce = exp(t(design.reduce %*% Beta.reduce))
    V.reduce = sweep(U.reduce, 2, k, FUN = "*")
    
    
    list(Beta.reduce = t(Beta.reduce), mean.reduce = U.reduce, dispersion = phi, 
         df.reduce = nbeta)
}


glmfit.OneGroup <- function(NanoStringData, m0, sigma, phi) {
    
    c = positiveFactor(NanoStringData)
    d = housekeepingFactor(NanoStringData)
    k = c * d
    lamda_i = negativeFactor(NanoStringData)
    Y = exprs(NanoStringData)
    Y_n = sweep(Y, 2, lamda_i, FUN = "-")
    Y_nph = sweep(Y_n, 2, k, FUN = "/")
    Y_nph[Y_nph <= 0] = 0.1
    
    
    
    n = ncol(Y)
    max.mean = max(rowMeans(Y_nph))
    
    get.mu <- function(dat) {
        
        y = dat[1:n]
        phi = dat[n + 1]
        
        obj = function(mu) {
            
            alpha = 1/phi
            tmp1 = 1/(1 + mu * k * phi)
            tmp2 = 1 - tmp1
            tmp2[tmp2==0] = 1e-08
            
            
            item1 = function(yy) {
                y_gi = yy[1]
                lamda_gi = yy[2]
                tmp2_gi = yy[3]
                t = c(0:y_gi)
                com = matrix(700, length(t), 1)
                
                
                tmp33.t = exp(rowMins(cbind(lgamma(t + alpha) + (y_gi - t) * log(lamda_gi) + 
                                                t * log(tmp2_gi) - lfactorial(t) - lfactorial(y_gi - t), com)))
                tmp33.tt = log(max(sum(tmp33.t), 1e-08))
            }
            
            tmp3 = apply(cbind(matrix(y, ncol = 1), matrix(lamda_i, ncol = 1), matrix(tmp2, 
                                                                                      ncol = 1)), 1, item1)
            
            
            -(sum(tmp3) - n * lgamma(alpha) + alpha * sum(log(tmp1)) - ((log(phi) - 
                                                                             m0)^2)/(2 * (sigma^2)) - log(sigma) - sum(lamda_i))
        }
        
        
        return(optimize(obj, interval = c(0.1, max.mean))$minimum)
    }
    
    mu = apply(cbind(matrix(Y, ncol = n), matrix(phi, ncol = 1)), 1, get.mu)
    
    Beta.reduce = log(mu)
    U.reduce = matrix(rep(mu, n), ncol = n)
    V.reduce = sweep(U.reduce, 2, k, FUN = "*")
    eta = log(phi)
    
    list(Beta.reduce = Beta.reduce, mean.reduce = U.reduce, dispersion = phi, df.reduce = 1)
} 
