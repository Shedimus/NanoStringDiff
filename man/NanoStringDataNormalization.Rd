\name{NanoStringDataNormalization}
\alias{NanoStringDataNormalization,NanoStringDataNormalization-method}
\alias{NanoStringDataNormalization}
\title{
   Normalize NanoStringData
}
\description{
   This function is used to get Normalized NanoString Data after adjusting for positive size factors, background noise and housekeeping size factors. Note that the normalized data values should only be used for data exploration / visualization purposes, e.g. drawing a heatmap. To perform differential expression analysis, we recommend users to follow the procedure described in the package vignette.
}
\usage{
   NanoStringDataNormalization(path=path, header=TRUE, designs)
}

\arguments{
  \item{path}{ the path of the file which the data are to be read from.}
  \item{header}{ a logical value indicating whether the file contains the names of the variables as its first line. If missing, the value is determined from the file format: header is set to TRUE if and only if the first row contains one fewer field than the number of columns.}
  \item{designs}{ a data frame in which the length of vector matches the column number of NanoStringData}
}


\author{
  hong wang <hong.wang@uky.edu>
  tingting zhai <tingting.zhai@uky.edu>
  chi wang <chi.wang@uky.edu>
}

\examples{
##path="/Users/NanoStringdiff-Rcode/Data/horbinski.csv"
##designs=data.frame(control=c(0,0,0,1,1,1))
##NanoStringDataNormalization(path=path, header=TRUE, designs)
}

\keyword{ normalization }



