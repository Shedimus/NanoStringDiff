## generate count matrix
endogenous = matrix(c(100,100,100,100,100,100,50,50,300,300,100,100),
                    3,4, byrow = TRUE)
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
pheno=pData(NanoStringData)
group=pheno$group
design.full=model.matrix(~0+group)

result = glm.LRT(NanoStringData, design.full, contrast = c(1, -1))


lfc1 = result$table[1,1]
lfc2 = result$table[2,1]
lfc3 = result$table[3,1]
expect_equal(lfc1, log2(1), tolerance=.01)
expect_equal(lfc2, log2(2), tolerance=.2)
expect_equal(lfc3, log2(3), tolerance=.2)



result = glm.LRT(NanoStringData, design.full, contrast = c(-1, 1))

lfc1 = result$table[1,1]
lfc2 = result$table[2,1]
lfc3 = result$table[3,1]
expect_equal(lfc1, -log2(1), tolerance=.01)
expect_equal(lfc2, -log2(2), tolerance=.2)
expect_equal(lfc3, -log2(3), tolerance=.2)
