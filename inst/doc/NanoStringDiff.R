### R code from vignette source 'NanoStringDiff.Rnw'

###################################################
### code chunk number 1: quick start (eval = FALSE)
###################################################
## library("Biobase")
## library("NanoStringDiff")
## data(NanoStringData)
## pheno=pData(NanoStringData)
## design.full=model.matrix(~0+pheno$group)
## NanoStringData=estNormalizationFactors(NanoStringData)
## result=glm.LRT(NanoStringData,design.full,contrast=c(1,-1))


###################################################
### code chunk number 2: directory (eval = FALSE)
###################################################
## directory <- "/path/to/your/files/"


###################################################
### code chunk number 3: GetDirectory
###################################################
directory <- system.file("extdata", package="NanoStringDiff", mustWork=TRUE)
path<-paste(directory,"Mori.csv",sep="/")


###################################################
### code chunk number 4: MakePhenotypeData
###################################################
designs=data.frame(group=c("Normal","Normal","Tumor","Tumor"))
designs


###################################################
### code chunk number 5: CsvInput
###################################################
library("NanoStringDiff")
NanoStringData=createNanoStringSetFromCsv(path,header=TRUE,designs)
NanoStringData


###################################################
### code chunk number 6: MatrixInput
###################################################
endogenous=matrix(rpois(100,50),25,4)
colnames(endogenous)=paste("Sample",1:4)
positive=matrix(rpois(24,c(128,32,8,2,0.5,0.125)*80),6,4)
colnames(positive)=paste("Sample",1:4)
negative=matrix(rpois(32,10),8,4)
colnames(negative)=paste("Sample",1:4)
housekeeping=matrix(rpois(12,100),3,4)
colnames(housekeeping)=paste("Sample",1:4)
designs=data.frame(group=c("Control","Control","Treatment","Treatment"),
                   gender=c("Male","Female","Female","Male"),
                   age=c(20,40,39,37))
NanoStringData1=createNanoStringSet(endogenous,positive,
                                 negative,housekeeping,designs)
NanoStringData1
pData(NanoStringData1)
head(exprs(NanoStringData1))


###################################################
### code chunk number 7: MakeDesignMatrix1
###################################################
pheno=pData(NanoStringData1)
group=pheno$group
design.full=model.matrix(~0+group)
design.full


###################################################
### code chunk number 8: MakeContrast
###################################################
contrast=c(-1,1)


###################################################
### code chunk number 9: GetControls
###################################################
NanoStringData1=estNormalizationFactors(NanoStringData1)
positiveFactor(NanoStringData1)
negativeFactor(NanoStringData1)
housekeepingFactor(NanoStringData1)


###################################################
### code chunk number 10: result1
###################################################
result=glm.LRT(NanoStringData1,design.full,contrast=contrast)
head(result$table)
str(result)


###################################################
### code chunk number 11: CreatePseudoData
###################################################
endogenous=matrix(rpois(300,50),25,12)
colnames(endogenous)=paste("Sample", 1:12)
colnames(endogenous)=paste("Sample",1:12)
positive=matrix(rpois(72,c(128,32,8,2,0.5,0.125)*80),6,12)
negative=matrix(rpois(96,10),8,12)
housekeeping=matrix(rpois(36,100),3,12)
designs=data.frame(group=c(rep("A",4),rep("B",4),rep("C",4)),
                   gender=rep(c("Male","Male","Female","Female"),3),
                   age=c(20,40,39,37,29,47,23,45,34,65,35,64))
NanoStringData2=createNanoStringSet(endogenous,positive,
                                 negative,housekeeping,designs)
NanoStringData2
pData(NanoStringData2)


###################################################
### code chunk number 12: MakeDesignMatrix2
###################################################
pheno=pData(NanoStringData2)
group=pheno$group
design.full=model.matrix(~0+group)
design.full


###################################################
### code chunk number 13: AssignSizeFactors
###################################################
NanoStringData2=estNormalizationFactors(NanoStringData2)


###################################################
### code chunk number 14: result2 (eval = FALSE)
###################################################
## result=glm.LRT(NanoStringData2,design.full,contrast=c(-1,1,0))


###################################################
### code chunk number 15: result3 (eval = FALSE)
###################################################
## result=glm.LRT(NanoStringData2,design.full,contrast=c(-1,0,1))


###################################################
### code chunk number 16: result4 (eval = FALSE)
###################################################
## result=glm.LRT(NanoStringData2,design.full,contrast=c(0,1,-1))


###################################################
### code chunk number 17: MakeDesignMatrix3
###################################################
design.full=model.matrix(~group)
design.full


###################################################
### code chunk number 18: result5 (eval = FALSE)
###################################################
## result=glm.LRT(NanoStringData2,design.full,Beta=2)


###################################################
### code chunk number 19: result6 (eval = FALSE)
###################################################
## result=glm.LRT(NanoStringData2,design.full,Beta=3)


###################################################
### code chunk number 20: MakeDesignMatrix4 (eval = FALSE)
###################################################
## design.full=model.matrix(~0+group)
## result=glm.LRT(NanoStringData2,design.full,contrast=c(1,-1/2,-1/2))


###################################################
### code chunk number 21: MakeDesignMatrix5
###################################################
design=pData(NanoStringData2)[,c("group","gender")]
group=design$group
gender=design$gender
design.full=model.matrix(~group+group:gender)
design.full


###################################################
### code chunk number 22: result7 (eval = FALSE)
###################################################
## result=glm.LRT(NanoStringData2,design.full,Beta=4)


###################################################
### code chunk number 23: result8 (eval = FALSE)
###################################################
## result=glm.LRT(NanoStringData2,design.full,Beta=2)


###################################################
### code chunk number 24: result9 (eval = FALSE)
###################################################
## result=glm.LRT(NanoStringData2,design.full,Beta=c(2,5))


###################################################
### code chunk number 25: MakeDesignMatrix6
###################################################
design.full=model.matrix(~group+gender+group:gender)


###################################################
### code chunk number 26: MakeDesignMatrix7
###################################################
design.full=model.matrix(~group*gender)
design.full


###################################################
### code chunk number 27: result10 (eval = FALSE)
###################################################
## result=glm.LRT(NanoStringData2,design.full,Beta=4:5)


###################################################
### code chunk number 28: sessionInfo
###################################################
toLatex(sessionInfo())


