\name{chipSmooth<-}


\alias{chipSmooth<-}

\title{
Setter method for \code{chipSmooth} slot in
\code{\link{parameterOptions}} object.
}
\description{
Setter method for \code{chipSmooth} slot in
\code{\link{parameterOptions}} object.}
\usage{
chipSmooth(object) <- value
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{parameterOptions}} object.
    }
    \item{value}{
    \code{value} is the positive numeric value to be assigned to the
    \code{chipSmooth} slot in \code{\link{parameterOptions}}
    Default value is 250 base pairs.
    }
}
\details{
When computing ChIP-seq like (\code{\link{computeChIPProfile}}) profile
from occupancy data (see \code{\link{computeOccupancy}}),
the profiles are smoothed using a window of a given size.
The default value is set at 250 base pairs.If \code{chipSmooth}
is set to 0 then the profile will not be smoothed.
We strongly encourage to set values when building a
\code{\link{parameterOptions}} object.
}
\value{
Returns a \code{\link{parameterOptions}} object with an updated
value for \code{chipSmooth} slot.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.

}
\author{
Patrick C.N Martin <pm16057@essex.ac.uk>
}

\examples{
# Building parameterOptions object
OPP <- parameterOptions()
# Setting new value for chipSd slot
chipSmooth(OPP) <-  250
}
