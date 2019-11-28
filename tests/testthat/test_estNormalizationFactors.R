## generate count matrix
endogenous = matrix(rpois(100,50),25,4)
sample.names = paste0("Sample",1:4)
colnames(endogenous) = sample.names
positive = matrix(c(128,32,8,2,0.5,0.125)*80,6,4)
colnames(positive) = paste("Sample",1:4)
negative = matrix(10,8,4)
colnames(negative) = sample.names
housekeeping = matrix(100,3,4)
colnames(housekeeping) = sample.names

## generate phenotype data
designs = data.frame(group=c("Control","Control","Treatment","Treatment"),
                   gender=c("Male","Female","Female","Male"),
                   age=c(20,40,39,37))
## input data to create a "NanoStringSet" object
NanoStringData = createNanoStringSet(endogenous,positive,
                                negative,housekeeping,designs)
NanoStringData = estNormalizationFactors(NanoStringData)




posFactor = matrix(1,1,4)
houFactor = matrix(1,1,4)

expect_true(all(positiveFactor(NanoStringData) == posFactor))
expect_true(all(housekeepingFactor(NanoStringData) == houFactor))








