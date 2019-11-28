\name{PlotsPositiveHousekeeping}
\alias{PlotsPositiveHousekeeping,PlotsPositiveHousekeeping-method}
\alias{PlotsPositiveHousekeeping}
\title{
   Plots of positive controls and housekeeping genes
}
\description{
   This function is used to pre-check the expressions of positive controls and housekeeping genes before data analysis. Linear regression plot of positive controls and variation analysis of housekeeping genes are available. The expressions of positive controls are supposed to be linearly related to the concentration of input sample materials, and the expressions of housekeeping genes are supposed to have relatively low variation. Nanostring recommends at least three housekeeping genes, but the more that are included, the more accurate the normalization will be. 
}
\usage{
   PlotsPositiveHousekeeping(path=path, header=TRUE)
}

\arguments{
  \item{path}{ the path of the file which the data are to be read from.}
  \item{header}{ a logical value indicating whether the file contains the names of the variables as its first line. If missing, the value is determined from the file format: header is set to TRUE if and only if the first row contains one fewer field than the number of columns.}
}


\author{
  hong wang <hong.wang@uky.edu>
  tingting zhai <tingting.zhai@uky.edu>
  chi wang <chi.wang@uky.edu>
}

\examples{
##path="/Users/NanoStringdiff-Rcode/Data/horbinski.csv"
##PlotsPositiveHousekeeping(path=path, header=TRUE)
}

\keyword{Plots, Positive Controls, Housekeeping Genes}



