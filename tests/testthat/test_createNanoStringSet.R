## generate count matrix
endogenous = matrix(rpois(100,50),25,4)
sample.names = paste0("Sample",1:4)
colnames(endogenous) = sample.names
positive = matrix(rpois(24,c(128,32,8,2,0.5,0.125)*30),6,4)
colnames(positive) = paste("Sample",1:4)
negative = matrix(rpois(32,10),8,4)
colnames(negative) = sample.names
housekeeping = matrix(rpois(12,100),3,4)
colnames(housekeeping) = sample.names

## generate phenotype data
designs = data.frame(group=c("Control","Control","Treatment","Treatment"),
                   gender=c("Male","Female","Female","Male"),
                   age=c(20,40,39,37))
## input data to create a "NanoStringSet" object
NanoStringData = createNanoStringSet(endogenous,positive,
                                 negative,housekeeping,designs)

expect_true(all(exprs(NanoStringData) == endogenous))
expect_true(all(positiveControl(NanoStringData) == positive))
expect_true(all(negativeControl(NanoStringData) == negative))
expect_true(all(housekeepingControl(NanoStringData) == housekeeping))
expect_true(all(colnames(exprs(NanoStringData)) == sample.names))








