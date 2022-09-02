\name{backgroundSignal}


\alias{backgroundSignal}


\title{
Accessor method for the \code{backgroundSignal} slot in a
\code{\link{parameterOptions}} object.
}
\description{
Extract or access the \code{backgroundSignal} slot in a
\code{\link{parameterOptions}} object.
}
\usage{
backgroundSignal(object)
}

\arguments{
    \item{object}{
    \code{object} is an \code{\link{parameterOptions}}
    }
}
\details{
Default Value: 0

When computing \code{\link{computeOccupancy}}, a ChIP-seq background
signal is used to scale Occupancy by considering both a
\code{backgroundSignal} and a \code{\link{maxSignal}}.
The \code{backgroundSignal}
is also used to nomalise occupancies against maxOccupancy.
The \code{backgroundSignal} usually comes from
experimental data and is provided by user. As a general rule,
if ChIP-seq data is available and will be used in
\code{\link{computeChIPProfile}}
, \code{\link{profileAccuracyEstimate}} or \code{\link{plotOccupancyProfile}},
it is advised to use the \code{backgroundSignal} from this data.
We strongly encourage to set values when building a
\code{\link{parameterOptions}} object.
}
\value{
Returns a \code{backgroundSignal} of a
\code{\link{parameterOptions}} object.
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
#Viewing single value in object
backgroundSignal(OPP)


}
