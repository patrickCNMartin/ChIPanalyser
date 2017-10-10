\name{PositionWeightMatrix}
\alias{PositionWeightMatrix}
\alias{PositionWeightMatrix-methods}

\title{
Accessor Method for the \code{PWM} slot in a
\code{\link{genomicProfileParameters}} object
}
\description{
Accessor Method for the \code{PWM} slot in a
\code{\link{genomicProfileParameters}} object
}
\usage{
PositionWeightMatrix(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfileParameters}}
    }
}
\details{
After creating a \code{\link{genomicProfileParameters}} object,
it is possible to access the Position Weight Matrix stored in this slot.
This slot should always contain something. This slot is either supplied by
user or directly computed from a Position Frequency Matrix when supplied.
}
\value{
Returns a matrix in the form of a Position Weight Matrix
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>
}

\examples{
#Loading data
data(ChIPanalyserData)
#Loading PFM files
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BCDSlx.pfm")
#Building genomicProfileParameters object
GPP<-genomicProfileParameters(PFM=PFM)
# Accessing Slot
PositionWeightMatrix(GPP)
}
