setClass("NanoStringSet", 
         contains = "ExpressionSet", 
         representation = representation(
             positiveFactor = "numeric", 
             negativeFactor = "numeric", 
             housekeepingFactor = "numeric", 
             positiveControl = "matrix", 
             negativeControl = "matrix", 
             housekeepingControl = "matrix"))
### normalization factor positiveFactor
setGeneric("positiveFactor", function(object) standardGeneric("positiveFactor"))
setGeneric("positiveFactor<-", function(object, value) standardGeneric("positiveFactor<-"))
setMethod("positiveFactor", signature(object = "NanoStringSet"), function(object) {
    object@positiveFactor
    
})
setReplaceMethod("positiveFactor", signature(object = "NanoStringSet", value = "numeric"), 
                 function(object, value) {
                     n = ncol(exprs(object))
                     if (length(value) != n) 
                         stop("wrong length for positive factor vector!")
                     object@positiveFactor <- value
                     object
                 })
### negativeFactor
setGeneric("negativeFactor", function(object) standardGeneric("negativeFactor"))
setGeneric("negativeFactor<-", function(object, value) standardGeneric("negativeFactor<-"))
setMethod("negativeFactor", signature(object = "NanoStringSet"), function(object) {
    object@negativeFactor
})
setReplaceMethod("negativeFactor", signature(object = "NanoStringSet", value = "numeric"), 
                 function(object, value) {
                     n = ncol(exprs(object))
                     if (length(value) != n) 
                         stop("wrong length for negative factor vector!")
                     object@negativeFactor <- value
                     object
                 })
### housekeepingFactor
setGeneric("housekeepingFactor", function(object) standardGeneric("housekeepingFactor"))
setGeneric("housekeepingFactor<-", function(object, value) standardGeneric("housekeepingFactor<-"))
setMethod("housekeepingFactor", signature(object = "NanoStringSet"), function(object) {
    object@housekeepingFactor
})
setReplaceMethod("housekeepingFactor", signature(object = "NanoStringSet", value = "numeric"), 
                 function(object, value) {
                     n = ncol(exprs(object))
                     if (length(value) != n) 
                         stop("wrong length for housekeeping factor vector!")
                     object@housekeepingFactor <- value
                     object
                 })
################## control genes positive control genes
setGeneric("positiveControl", function(object) standardGeneric("positiveControl"))
setGeneric("positiveControl<-", function(object, value) standardGeneric("positiveControl<-"))
setMethod("positiveControl", signature(object = "NanoStringSet"), function(object) {
    object@positiveControl
})
setReplaceMethod("positiveControl", signature(object = "NanoStringSet", value = "matrix"), 
                 function(object, value) {
                     object@positiveControl <- value
                     object
                 })
### negative control genes
setGeneric("negativeControl", function(object) standardGeneric("negativeControl"))
setGeneric("negativeControl<-", function(object, value) standardGeneric("negativeControl<-"))
setMethod("negativeControl", signature(object = "NanoStringSet"), function(object) {
    object@negativeControl
})
setReplaceMethod("negativeControl", signature(object = "NanoStringSet", value = "matrix"), 
                 function(object, value) {
                     object@negativeControl <- value
                     object
                 })
### housekeeping control genes
setGeneric("housekeepingControl", function(object) standardGeneric("housekeepingControl"))
setGeneric("housekeepingControl<-", function(object, value) standardGeneric("housekeepingControl<-"))
setMethod("housekeepingControl", signature(object = "NanoStringSet"), function(object) {
    object@housekeepingControl
})
setReplaceMethod("housekeepingControl", signature(object = "NanoStringSet", value = "matrix"), 
                 function(object, value) {
                     object@housekeepingControl <- value
                     object
                 })
################# read data#####
createNanoStringSet <- function(endogenous, positiveControl, 
                             negativeControl, housekeepingControl, designs) {
    ## work on input design matrix. Carefully about the input design.
    if (is(designs, "vector")) {
        designs <- as.data.frame(designs)
        colnames(designs) = "designs"
    } else if (is(designs, "matrix")) {
        ## multiple factor design still single factor
        if (ncol(designs) == 1) {
            designs <- as.data.frame(designs)
            colnames(designs) = "designs"
        } else {
            ## input is a data frame
            designs <- as.data.frame(designs)
        }
    }
    
    rownames(designs) = colnames(endogenous)
    
    
    if (is(designs, "data.frame") || is(designs, "AnnotatedDataFrame")) {
        stopifnot(nrow(designs) == ncol(endogenous))
        designs <- as(designs, "AnnotatedDataFrame")
    }
    res <- new("NanoStringSet", exprs = endogenous, phenoData = designs, 
               positiveControl = positiveControl, 
               negativeControl = negativeControl, 
               housekeepingControl = housekeepingControl)
    
    res
    
}
createNanoStringSetFromCsv <- function(path, header = TRUE, designs) {
    
    ## read data from csv file
    
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
    spikein = counts[spikein.id, ]
    endogenous = counts[-c(pos.id, neg.id, house.id, spikein.id), ]
    
    
    cat(paste(" There are", ncol(counts), "samples imported; \n There are ", 
              nrow(counts), "genes imported with:"))
    
    
    print(table(code.class))
    
    
    
    
    ## work on input design matrix. Carefully about the input design.
    if (is(designs, "vector")) {
        designs <- as.data.frame(designs)
        colnames(designs) = "designs"
    } else if (is(designs, "matrix")) {
        ## multiple factor design still single factor
        if (ncol(designs) == 1) {
            designs <- as.data.frame(designs)
            colnames(designs) = "designs"
        } else {
            ## input is a data frame
            designs <- as.data.frame(designs)
        }
    }
    
    
    rownames(designs) = colnames(counts)
    
    
    if (is(designs, "data.frame") || is(designs, "AnnotatedDataFrame")) {
        stopifnot(nrow(designs) == ncol(counts))
        designs <- as(designs, "AnnotatedDataFrame")
    }
    
    
    
    rownames(designs) = colnames(counts)
    
    
    
    res <- new("NanoStringSet", exprs = endogenous, phenoData = designs, 
               positiveControl = positive, 
               negativeControl = negative, 
               housekeepingControl = housekeeping)
    
    res
    
} 
