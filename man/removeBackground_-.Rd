\name{removeBackground<-}
\alias{removeBackground<-}
\alias{removeBackground<--methods}


\title{
Setter Method for the \code{removeBackground} slot in a
\code{occupancyProfileParameters} object
}
\description{
Setter Method for the \code{removeBackground} slot in a
\code{occupancyProfileParameters} object
}
\usage{
removeBackground(object) <-value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{
    \code{object} is  an \code{\link{occupancyProfileParameters}} object
    }
    \item{value}{
    \code{value} is positive numerical value to be assigned to the
    \code{removeBackground} slot in a \code{\link{occupancyProfileParameters}}
    object. Default is set a 0.
    }
}
\details{
A numeric value describing a threshold at which Occupancy signals must be
removed (Default is set at 0). The removal of Occupancy signals will occur
when computing \code{\link{computeOccupancy}}
(see \code{\link{computeOccupancy}} function)
}
\value{
Returns an \code{\link{occupancyProfileParameters}} object  with an updated
value for the \code{removeBackground} slot
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>
}

\examples{
#Building occupancyProfileParameters object
OPP <- occupancyProfileParameters()
#Setting new Value for removeBackground
removeBackground(OPP) <- 0.1
}
