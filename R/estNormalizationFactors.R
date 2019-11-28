estNormalizationFactors <- function(NanoStringData) {
    
    if (!is(NanoStringData, "NanoStringSet")) 
        stop("Input must be an object of NanoStringSet class!")
    
    
    
    pos = positiveControl(NanoStringData)
    pos[pos <= 0] = 1
    neg = negativeControl(NanoStringData)
    x = c(128, 32, 8, 2, 0.5, 0.125, rep(0, nrow(neg)))
    Y = rbind(pos, neg)
    n = ncol(Y)
    lamda_i = rep(0, n)
    c = rep(0, n)
    for (i in 1:n) {
        
        model = glm(Y[, i] ~ x, family = poisson(link = identity))
        lamda_i[i] = round(model$coeff[1])
        
        c[i] = model$coeff[2]
    }
    
    coeff = c  # real coefficients
    c = c/mean(c)  # size factors
    
    names(lamda_i) = colnames(pos)
    names(c) = colnames(pos)
    lamda_i[lamda_i == 0] = 1
    
    
    
    ## estimate housekeeping factors
    house1 = sweep(housekeepingControl(NanoStringData), 2, lamda_i, FUN = "-")
    house1[house1 <= 0] = 1
    house2 = sweep(house1, 2, c, FUN = "/")
    temp2 = sweep(house2, 1, rowMeans(house2), FUN = "/")
    d = apply(temp2, 2, median)  #housekeeping control factors
    
    positiveFactor(NanoStringData) = c
    negativeFactor(NanoStringData) = lamda_i
    housekeepingFactor(NanoStringData) = d
    NanoStringData
} 
