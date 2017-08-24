\name{chipSmooth}
\alias{chipSmooth}
\alias{chipSmooth-methods}

\title{
Accessor methods for \code{chipSmooth} slot in a
\code{\link{occupancyProfileParameters}} object.
}
\description{
Access or Extract \code{chipSmooth} slot in a
\code{\link{occupancyProfileParameters}} object.
}
\usage{
chipSmooth(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{occupancyProfileParameters}} object.
    }
}
\details{
When computing ChIP-seq like (\code{\link{computeChipProfile}}) profile
from occupancy data (see \code{\link{computeOccupancy}}),
the profiles are smoothed using a window of a given size.
The default value is set at 250 base pairs. If \code{chipSmooth}
is set to 0 then the profile will not be smoothed.
We strongly encourage to set values when building a
\code{\link{occupancyProfileParameters}} object.
}
\value{
Returns the \code{chipSmooth} slot in an
\code{\link{occupancyProfileParameters}} object.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}
\author{
Patrick C.N. Martin <pm16057@essex.ac.uk>
}
\examples{
# Building occupancyProfileParameters object
OPP <- occupancyProfileParameters()
# Accessing chipSd slot
chipSmooth(OPP)
}
