\name{positiveControl}
\alias{positiveControl,NanoStringSet-method}
\alias{positiveControl<-,NanoStringSet,matrix-method}
\alias{positiveControl}
\alias{positiveControl<-}
\title{
Acessor functions for the 'positiveControl' slot in a
  NanoStringSet object.
}
\description{
  nCounter Analyzer has positive spike-in RNA hybridization controls 
  for each sample which can be used to estimate the overall efficiency 
  of hybridization and recovery for each sample.

}
\usage{
\S4method{positiveControl}{NanoStringSet}(object)
\S4method{positiveControl}{NanoStringSet,matrix}(object) <- value
}

\arguments{

 \item{object}{A NanoStringSet object.}
 \item{value}{A matrix with six positive control genes.}
}

\value{A matrix contain positive control genes}

\details{
  Positive control genes are provided by NanoString nCounter technology.
  For each sample, nCounter provide six positive controls corresponding 
  to six different concentrations in the 30 ul hybridzation:
  128fM, 32fM, 8fM, 2fM, 0.5fM, and 0.125fM. Six positive control genes must 
  be order by concentrations from high to low.
}
\author{
   Hong Wang <hong.wang@uky.edu>
   chi wang <chi.wang@uky.edu>
}

\seealso{
   positiveFactor
}
\examples{
data(NanoStringData)
## obtain positive control genes
positiveControl(NanoStringData)

## assign a matrix
n=ncol(exprs(NanoStringData))
x=matrix(c(128,32,8,2,0.5,0.125)*80,ncol=1)
positive=matrix(rpois(6*n,x),ncol=n)
positiveControl(NanoStringData)=positive
}

