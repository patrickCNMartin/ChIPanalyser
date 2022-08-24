\name{backgroundSignal<-}


\alias{backgroundSignal<-}
\title{
Setter method for \code{backgroundSignal} slot in a
\code{\link{parameterOptions}}
}
\description{
Setter method for \code{backgroundSignal} slot in a
\code{\link{parameterOptions}}

}
\usage{
backgroundSignal(object)<-value
}

\arguments{
    \item{object}{
    \code{object} is an \code{\link{parameterOptions}} object.
    }
    \item{value}{
    \code{value} is the value to be assigned to the \code{backgroundSignal}
    slot in
    \code{\link{parameterOptions}}.
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
\code{\link{computeChIPProfile}}
, \code{\link{profileAccuracyEstimate}} or
\code{\link{plotOccupancyProfile}}, it is advised to use
the \code{backgroundSignal} from this data.
We strongly encourage to set values when building a
\code{\link{parameterOptions}} object.
}
\value{
Returns a \code{\link{parameterOptions}} object with a new
\code{value} assigned to the \code{backgroundSignal} slot.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
Patrick C.N. Martin and Nicolae Radu Zabe (2020) Dissecting the binding mechanisms of transcription factors to DNA using a statistical thermodynamics framework. CSBJ, 18, 3590-3605.
}
\author{
Patrick C.N. Martin <pcnmartin@gmail.com>
}

\examples{
# Building occupancyProfileParameters object
OPP <- parameterOptions()
# Setting new value for backgroundSignal
backgroundSignal(OPP) <- 0.2
# Viewing whole object with new updated value
OPP
#Viewing single value in object
backgroundSignal(OPP)

}
