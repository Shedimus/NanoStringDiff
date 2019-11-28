\name{negativeControl}
\alias{negativeControl,NanoStringSet-method}
\alias{negativeControl<-,NanoStringSet,matrix-method}
\alias{negativeControl}
\alias{negativeControl<-}
\title{
Acessor functions for the 'negativeControl' slot in a
  NanoStringSet object.
}
\description{
Negative control genes are provided by nCounter Analyzer which can be used to 
estimate background noise for each sample.
}
\usage{
\S4method{negativeControl}{NanoStringSet}(object)
\S4method{negativeControl}{NanoStringSet,matrix}(object) <- value
}

\arguments{

 \item{object}{A NanoStringSet object.}
 \item{value}{A matrix with negative control genes.}
}

\value{A matrix contain negative control genes}

\details{
  Each code set in the nCounter Analyzer includes several negatives 
  control genes for which no tranCounterript is expected to be present.
  We use these spike-in negative control genes to estimate
  background noise for each sample.
}
\author{
   Hong Wang <hong.wang@uky.edu>
   chi wang <chi.wang@uky.edu>
}

\seealso{
   negativeFactor
}
\examples{
data(NanoStringData)
## obtain negative control genes
negativeControl(NanoStringData)

## assign a matrix
n=ncol(exprs(NanoStringData))
r=nrow(negativeControl(NanoStringData))
negative=matrix(rpois(r*n,10),ncol=n)
negativeControl(NanoStringData)=negative
}

