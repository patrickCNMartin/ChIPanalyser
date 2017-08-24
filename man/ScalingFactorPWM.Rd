\name{ScalingFactorPWM}
\alias{ScalingFactorPWM}
\alias{ScalingFactorPWM-methods}


\title{
Accessor method for the \code{lambda} slot in  a
\code{\link{genomicProfileParameters}} object
}
\description{
Accessor method for the \code{lambda} slot in  a
\code{\link{genomicProfileParameters}} object
}
\usage{
ScalingFactorPWM(object)
}
\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfileParameters}} object
    }
}
\details{
The model used in \code{ChIPanalyser} relies on a scaling factor
(scaling Transcription Factor binding affinity) for the PWM Score
(or binding energy). It is also referred to as \code{lambda}.
In the case a vector of values is provided, each value will be either
computed individually ( as it will be the case in
\code{\link{computeGenomeWidePWMScore}}) or in combination with
\code{\link{boundMolecules}} (as it will be the case in
\code{\link{computeOccupancy}}).
}
\value{
Returns the value or vector assigned to the \code{lambda}
slot in a \code{\link{genomicProfileParameters}} object
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>
}
\examples{
# Loading data
data(ChIPanalyserData)
#Loading PFM files
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BCDSlx.pfm")
#Building data objects
GPP <- genomicProfileParameters(PFM=PFM, ScalingFactorPWM=c(0.5,1,1.5))
#Accessing value for ScalingFactorPWM
ScalingFactorPWM(GPP)
}
