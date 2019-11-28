PlotsPositiveHousekeeping<- function(path=path, header=TRUE){
  # ------------------------------------ #
  # path: points out the directory in which your csv.file is located #
  # header: a logical value(TRUE or FALSE) indicating whether the file contains the names of the variables as its first line. #
  # designs: data frame for pheno type data storage #
  # ------------------------------------ #
  data = read.table(path, header = header, sep = ",")
  
  if (is.null(data)) {
      stop("There is no counts data")
  }
  
  selectcol = !(names(data) %in% c("Code.Class", "Name", "Accession"))
  
  ## remove NA from data set
  counts = data[,selectcol]
  counts = as.matrix(counts)
  id = which(is.na(rowSums(counts)))
  if (length(id) > 0) {
      data = data[-id, ]
  }
  
  code.class = tolower(data[, c("Code.Class")])
  name = data[, c("Name")]
  accession = data[, c("Accession")]
  
  counts = data[,selectcol]
  counts = as.matrix(counts)
  rownames(counts) = name
  
  pos.id = grep("positive", code.class, fixed = TRUE)
  neg.id = grep("negative", code.class, fixed = TRUE)
  house.id = grep("housekeeping", code.class, fixed = TRUE)
  spikein.id = grep("spikein", code.class, fixed = TRUE)
  
  positive = counts[pos.id, ]
  negative = counts[neg.id, ]
  housekeeping = counts[house.id, ]
  house.adjust = counts[house.id, ]
  
  mean.neg<-colMeans(negative)
  
  pos = positive
  pos[pos <= 0] = 1
  neg = negative
  x = c(128, 32, 8, 2, 0.5, 0.125, rep(0, nrow(neg)))
  Y = rbind(pos, neg)
  n = ncol(Y)
  c = rep(0, n)
  for (i in 1:n) {
      
      model = glm(Y[, i] ~ x, family = poisson(link = identity))
      c[i] = model$coeff[2]
  }
  
  coeff = c  # real coefficients
  c = c/mean(c)  # size factors
  
  for (i in 1:nrow(housekeeping)){
      housekeeping[i,]<-(housekeeping[i,]-mean.neg)/c^(1/2)
  }
  
  nhouse<-nrow(housekeeping)
  varhouse<-rowVars(housekeeping)^(1/2)/rowMeans(housekeeping)
  con=c(128,32,8,2,0.5,0.125)
  par(mfrow=c(1,2),mar=c(5,6,2,1)+0.1, oma=c(2,2,2,2)+0.1)
  plot(con,rowMeans(positive),xlab="Concentration",ylab="Positive Counts",
  col = "blue")
  abline(lm(rowMeans(positive) ~ con))
  plot(c(1:nhouse),varhouse,xaxt='n',xlab="",ylab="Coefficient of Variation",
  col = "blue", xlim =c(1,nhouse+1), ylim =c(min(varhouse)*0.9,max(varhouse)*1.1))
  title(xlab="Housekeeping Genes", line=1.5)
  textxy(c(1:nhouse), varhouse, labs=rownames(housekeeping), cex = 1, m = c(0, 0), offset = 0.65)
} # end function #
