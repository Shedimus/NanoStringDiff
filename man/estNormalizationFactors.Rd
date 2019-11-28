\name{estNormalizationFactors}
\alias{estNormalizationFactors,estNormalizationFactors-method}
\alias{estNormalizationFactors}
\title{
   estimate normalization factors,include positive size factors, background 
   noise, housekeeping size factors.
}
\description{
   This function estimates positive size factors, background noise and 
   housekeeping size factors for the input "NanoStringSet" object and 
   return the same object with positiveFactor, negativeFactor and 
   housekeepingFactor slots filled or replaced.
}
\usage{
   estNormalizationFactors(NanoStringData)
}

\arguments{
  \item{NanoStringData}{ An object of "NanoStringSet" class.}
}

\value{
  The same "NanoStringSet" object with positiveFactor, negativeFactor and 
  housekeepingFactor field filled or replaced. 
}

\author{
  hong wang <hong.wang@uky.edu>
  chi wang <chi.wang@uky.edu>
}

\examples{
data(NanoStringData)
NanoStringData=estNormalizationFactors(NanoStringData)
pf=positiveFactor(NanoStringData)
nf=negativeFactor(NanoStringData)
hf=housekeepingFactor(NanoStringData)
}

\keyword{ normalization }



