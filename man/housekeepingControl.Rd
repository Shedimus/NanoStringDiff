\name{housekeepingControl}
\alias{housekeepingControl,NanoStringSet-method}
\alias{housekeepingControl<-,NanoStringSet,matrix-method}
\alias{housekeepingControl}
\alias{housekeepingControl<-}
\title{
Acessor functions for the 'housekeepingControl' slot in a
  NanoStringSet object.
}
\description{
 user-defined housekeeping control genes can be used to estimate housekeeping 
 factors to adjust variation caused by different sample input.
}
\usage{
\S4method{housekeepingControl}{NanoStringSet}(object)
\S4method{housekeepingControl}{NanoStringSet,matrix}(object) <- value
}

\arguments{

 \item{object}{A NanoStringSet object.}
 \item{value}{A matrix with housekeeping control genes.}
}

\value{A matrix contain housekeeping control genes}

\details{
  NanoString nCounter analyzer also contains probes for a set of 
  species-specific mRNA housekeeping(reference) genes that are not 
  spike-in the system. Nanostring recommends at least three housekeeping genes, 
  but the more that are included, the more accurate the normalization will be. 
  Housekeeping control genes are expected consistent in their expression levels.
}
\author{
   Hong Wang <hong.wang@uky.edu>
   chi wang <chi.wang@uky.edu>
}

\seealso{
   housekeepingFactor
}
\examples{
data(NanoStringData)
## obtain housekeeping control genes
housekeepingControl(NanoStringData)

## assign a matrix
n=ncol(exprs(NanoStringData))
r=nrow(housekeepingControl(NanoStringData))
housekeeping=matrix(rpois(r*n,1000),ncol=n)
housekeepingControl(NanoStringData)=housekeeping
}


