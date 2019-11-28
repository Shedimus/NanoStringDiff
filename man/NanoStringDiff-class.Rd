\name{NanoStringSet-class}
\Rdversion{1.1}
\docType{class}
\alias{NanoStringSet-class}
\alias{NanoStringSet}
% Constructor-like function:
\alias{createNanoStringSet}
\alias{createNanoStringSetFromCsv}
\title{NanoStringSet object and constructors}
\description{
 The \code{NanoStringSet} is a s4 class used to store data from NanoString 
 nCounter analyzer. This class a subclass of \code{ExpressionSet}, with six 
 more slots: positiveControl, negativeControl, housekeepingControl, 
 positiveFactor, negativeFactor and housekeepingFactor.
 The constructor functions \code{createNanoStringSet} and 
 \code{createNanoStringSetFromCsv} create a NanoStringSet object from two types 
 of input: seperate matrix or csv files. See the vignette for examples of 
 contruction from these two input types.
 }

\usage{
createNanoStringSet(endogenous,positiveControl,negativeControl,
                 housekeepingControl,designs)

createNanoStringSetFromCsv(path, header=TRUE, designs)
}

\arguments{
  \item{endogenous}{for matrix input: a matrix of non-negative integers 
                    of endogenes}
  
  \item{positiveControl}{for matrix input: a matrix of non-negative integers of 
        positive control genes. There must have 6 positive control genes order 
        by concentrations form high to low}
  
  \item{negativeControl}{for matrix input: a matrix of non-negative integers 
       of negative control genes}
        
  \item{housekeepingControl}{for matrix input: a matrix of non-negative 
       integers of housekeeping control genes}
        
  \item{designs}{for data.frame input: phenotype data for NanoString nCounter 
       data with at least one column. Each row is one sample, that is the 
       number of rows must equal number of samples or replicates in the data.}
  
  \item{path}{path to the csv file.}
  
  \item{header}{a logical value indicating whether the file contains the 
  names of the variables as its first line. The default value is TRUE.}
  
}

\value{
     A NanoStringSet object. 
}


\section{Methods}{
  \describe{
    \item{positiveControl, positiveControl<- :}{Access and set positive 
    control genes.}
    \item{negativeControl, negativeControl<- :}{Access and set negative 
    control genes.}
    \item{housekeepingControl, housekeepingControl<- :}{Access and set 
    housekeeping 
          control genes.}
    \item{positiveFactor, positiveFactor<- :}{Access and set positive factors.}
    \item{negativeFactor, negativeFactor<- :}{Access and set negative factors.}
    \item{housekeepingFactor, housekeepingFactor<- :}{Access and set 
    housekeeping factors.}
  }
}

\author{
    hong wang <hong.wang@uky.edu>
    chi wang <chi.wang@uky.edu>
}

\seealso{
  positiveControl, negativeControl, housekeepingControl,
  positiveFactor,  negativeFactor,  housekeepingFactor
}

\examples{
endogenous=matrix(rpois(100,50),25,4)
positive=matrix(rpois(24,c(128,32,8,2,0.5,0.125)*80),6,4)
negative=matrix(rpois(32,10),8,4)
housekeeping=matrix(rpois(12,100),3,4)
designs=data.frame(group=c(0,0,1,1),gender=c("male","female","female","male"),
age=c(20,40,39,37))
NanoStringData=createNanoStringSet(endogenous,positive,negative,
                                housekeeping,designs)
NanoStringData
pData(NanoStringData)
positiveControl(NanoStringData)
head(exprs(NanoStringData))

}

\keyword{ classes }

