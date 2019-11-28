\name{positiveFactor}
\alias{positiveFactor,NanoStringSet-method}
\alias{positiveFactor<-,NanoStringSet,numeric-method}
\alias{positiveFactor}
\alias{positiveFactor<-}
\title{
Acessor functions for the 'positiveFactor' slot in a
  NanoStringSet object.
}
\description{
  Positive size factors can be used to adjust all platform associated 
  sources of variation.

}
\usage{
\S4method{positiveFactor}{NanoStringSet}(object)
\S4method{positiveFactor}{NanoStringSet,numeric}(object) <- value
}

\arguments{

 \item{object}{A NanoStringSet object.}
 \item{value}{A vector of positive size factors.}
}

\value{A vector contain positive size factors}

\details{
 The observed counts including negative control genes and housekeeping control
 genes might be effect by some experimental factors like hybridization and 
 binding efficiency. In order to get the true rate of gene expression, 
 these variations must be normalized. Positive size factors can normalize 
 this kind of variation.
}
\author{
   Hong Wang <hong.wang@uky.edu>
   chi wang <chi.wang@uky.edu>
}

\seealso{
   positiveControl
}
\examples{
data(NanoStringData)
## obtain positive factors
positiveFactor(NanoStringData)

## assign a vector
n=ncol(exprs(NanoStringData))
positiveFactor(NanoStringData)=rep(1,n)
}

