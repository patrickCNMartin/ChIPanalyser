\name{backgroundSignal<-}
\alias{backgroundSignal<-}
\alias{backgroundSignal<--methods}

\title{
Setter method for \code{backgroundSignal} slot in a
\code{\link{occupancyProfileParameters}}
}
\description{
Setter method for \code{backgroundSignal} slot in a
\code{\link{occupancyProfileParameters}}

}
\usage{
backgroundSignal(object)<-value
}

\arguments{
    \item{object}{
    \code{object} is an \code{\link{occupancyProfileParameters}} object.
    }
    \item{value}{
    \code{value} is the value to be assigned to the \code{backgroundSignal}
    slot in
    \code{\link{occupancyProfileParameters}}.
    \code{backgroundSignal} should be a positive value.
    Defualt value is 0.
    }
}
\details{
Defualt value: 0.
When computing \code{\link{computeOccupancy}}, a ChIP-seq background
signal is used to
scale Occupancy by considering both a \code{backgroundSignal} and a
\code{\link{maxSignal}}. The \code{backgroundSignal}
is also used to nomalise occupancies to maxOccupancy.
The \code{backgroundSignal} usually comes from
experimental data and is provided by user. As a general rule,
if ChIP-seq data is available and will be used in
\code{\link{computeChipProfile}}
, \code{\link{profileAccuracyEstimate}} or
\code{\link{plotOccupancyProfile}}, it is advised to use
the \code{backgroundSignal} from this data.
We strongly encourage to set values when building a
\code{\link{occupancyProfileParameters}} object.
}
\value{
Returns a \code{\link{occupancyProfileParameters}} object with a new
\code{value} assigned to the \code{averageExpPWMScore} slot.
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
# Setting new value for backgroundSignal
backgroundSignal(OPP) <- 0.2
# Viewing whole object with new updated value
OPP
#Viewing single value in object
backgroundSignal(OPP)

}
