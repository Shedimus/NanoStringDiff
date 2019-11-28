\name{housekeepingFactor}
\alias{housekeepingFactor,NanoStringSet-method}
\alias{housekeepingFactor<-,NanoStringSet,numeric-method}
\alias{housekeepingFactor}
\alias{housekeepingFactor<-}
\title{
Acessor functions for the 'housekeepingFactor' slot in a
  NanoStringSet object.
}
\description{
  Housekeeping size factors can be used to adjust the variance 
  caused by different sample input.

}
\usage{
\S4method{housekeepingFactor}{NanoStringSet}(object)
\S4method{housekeepingFactor}{NanoStringSet,numeric}(object) <- value
}

\arguments{

 \item{object}{A NanoStringSet object.}
 \item{value}{A vector of housekeeping size factors.}
}

\value{A vector contain housekeeping factors}

\details{
 Housekeeping gene normalization corrects for different in sample input 
 between assays,since reference genes are suppose have same expression 
 rate between samples.So the read counts from housekeeping genes, after 
 subtracting background noise and adjusting by positive size factors, 
 that are not expected to vary between samples. If there exist difference, 
 which should be caused by sample input variation.
}
\author{
   Hong Wang <hong.wang@uky.edu>
   chi wang <chi.wang@uky.edu>
}

\seealso{
   housekeepingControl
}
\examples{
data(NanoStringData)
## obtain housekeeping factors
housekeepingFactor(NanoStringData)

## assign a vector
n=ncol(exprs(NanoStringData))
housekeepingFactor(NanoStringData)=rep(1,n)
}

