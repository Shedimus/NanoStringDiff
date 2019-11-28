\name{negativeFactor}
\alias{negativeFactor,NanoStringSet-method}
\alias{negativeFactor<-,NanoStringSet,numeric-method}
\alias{negativeFactor}
\alias{negativeFactor<-}
\title{
Acessor functions for the 'negativeFactor' slot in a
  NanoStringSet object.
}
\description{
  Negative size factors can be used to adjust background niose for each sample.

}
\usage{
\S4method{negativeFactor}{NanoStringSet}(object)
\S4method{negativeFactor}{NanoStringSet,numeric}(object) <- value
}

\arguments{

 \item{object}{A NanoStringSet object.}
 \item{value}{A vector of background noise.}
}

\value{A vector contain background noise}

\details{
 Accurate estimation of system background is essential for DE detection 
 analysis. Each code set in the nCounter Analyzer includes several negatives 
 control genes for which no tranCounterript is expected to be present.
 We use these spike-in negative control genes to estimate background 
 noise for each sample
}
\author{
   Hong Wang <hong.wang@uky.edu>
   chi wang <chi.wang@uky.edu>
}

\seealso{
   negativeControl
}
\examples{
data(NanoStringData)
## obtain negative factors
negativeFactor(NanoStringData)

## assign a vector
n=ncol(exprs(NanoStringData))
lamda=rpois(n,10)
negativeFactor(NanoStringData)=lamda
}

