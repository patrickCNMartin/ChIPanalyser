\name{removeBackground}
\alias{removeBackground}
\alias{removeBackground-methods}

\title{
Accessor Method for the \code{removeBackground} slot in a
\code{\link{occupancyProfileParameters}} object
}
\description{
Accessor Method for the \code{removeBackground} slot in a
\code{\link{occupancyProfileParameters}} object
}
\usage{
removeBackground(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{occupancyProfileParameters}} object
    }
}
\details{
A numeric value describing a threshold at which Occupancy signals must be
removed (Default is set at 0). The removal of Occupancy signals will occur
when computing \code{\link{computeOccupancy}}
(see \code{\link{computeOccupancy}} function)
}
\value{
Returns the value assigned to the \code{removeBackground} slot in a
\code{\link{occupancyProfileParameters}} object
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
#Accessing  Value for removeBackground
removeBackground(OPP)

}
