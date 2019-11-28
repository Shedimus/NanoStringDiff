\name{glm.LRT}
\alias{glm.LRT}
\alias{glm.LRT,NanoStringSet-method}

\title{
 perform gene-wise likelihood ratio test for NanoString Data
}
\description{

  The method considers a generalized linear model of the negative 
  binomial family to characterize count data and allows for multi-factor 
  design. The method propose an empirical Bayes shrinkage approach to 
  estimate the dispersion parameter and use likelihood ratio test to obtain
  p-value.
}
\usage{
   glm.LRT(NanoStringData,design.full,Beta=ncol(design.full), contrast=NULL)
}

\arguments{
  \item{NanoStringData}{An object of "NanoStringSet" class.}
  
  \item{design.full}{numeric matrix giving the design matrix for the 
                     generalized linear models under full model. must be of 
                     full column rank.}
  
  \item{Beta}{integer or character vector indicating which coefficients of the 
              linear model are to be tested equal to zero. Values must be 
              columns or column names of design. Defaults to the last 
              coefficient. Ignored if contrast is specified.}
  
  \item{contrast}{numeric vector or matrix specifying one or more contrasts 
              of the linear model coefficients to be tested equal to zero.}
}

\value{
  A list
 \item{table }{A data frame with each row corresponding to a gene. Rows are 
               sorted according to likelihood ratio test statistics. 
               The columns are:
               logFC: log fold change between two groups.
               lr: likelihood ratio test statictics.
               pvalue: p-value.
               qvalue: adjust p-value using the procedure of Benjamini 
               and Hochberg.}
               
 \item{dispersion }{a vertor of dispersion}
 
 \item{log.dispersion }{a vector of log dispersion: 
                        log.dispersion=log(dispersion)}
 
 \item{design.full}{numeric matrix giving the design matrix under full 
                    generalizedlinear model.}
                    
 \item{design.reduce}{numeric matrix giving the design matrix under reduced 
                      generalizedlinear model.}
 
 \item{Beta.full}{coefficients under full model.}
 
 \item{mean.full}{mean value under full model.}
 
 \item{Beta.reduce}{coefficients under reduced model.}
 
 \item{mean.reduce}{mean value under reduced model.}
 
 \item{m0}{hyper-parameter: mean value of the prior distribution of log 
           dispersion}
 
 \item{sigma}{hyper-parameter: standard deviation of the prior distribution of
              log dispersion}
}

\author{
   hong wang<hong.wang@uky.edu>
   chi wang <chi.wang@uky.edu>
}


\examples{
data(NanoStringData)
NanoStringData=estNormalizationFactors(NanoStringData)
group=pData(NanoStringData)
design.full=model.matrix(~0+factor(group$group))
contrast=c(1,-1)
result=glm.LRT(NanoStringData,design.full,
               Beta=ncol(design.full),contrast=contrast)
head(result$table)
}


\keyword{ models }
