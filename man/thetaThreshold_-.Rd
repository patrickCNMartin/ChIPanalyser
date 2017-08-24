\name{thetaThreshold<-}
\alias{thetaThreshold<-}
\alias{thetaThreshold<--methods}


\title{
Setter method the \code{thetaThreshold} slot in a
\code{\link{occupancyProfileParameters}} object
}
\description{
Setter method the \code{thetaThreshold} slot in a
\code{\link{occupancyProfileParameters}} object
}
\usage{
thetaThreshold(object) <- value
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{occupancyProfileParameters}}
    }
    \item{value}{
    \code{value} is a positive numeric value to be assigned to the
    \code{thetaThreshold} slot. Default is set at 0.1
    }
}
\details{
In order to measure the accuracy of predicted models agaisnt real
ChIP-seq data, \code{ChIPanalyser} uses three "measures": correlation,
Mean Squared Error and Theta. Theta is an in house metric that can be
described as a modified ratio of correlation over Mean Squarred Error.
The higher the value of theta , the better the fit becomes.
The theta threshold serves the purpose of avoiding extreme values and thus
clarifying optimal parameter heat maps
( see \code{\link{plotOptimalHeatMaps}}). Default is set at 0.1
}
\value{
Returns an \code{\link{occupancyProfileParameters}} object with an updated
value for \code{thetaThreshold}
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>
}

\examples{
# Building occupancyProfileParameters object
OPP <- occupancyProfileParameters()
#Setting new Value for thetaThreshold
thetaThreshold(OPP) <- 0.1
}
