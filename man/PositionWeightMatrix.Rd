\name{PositionWeightMatrix}


\alias{PositionWeightMatrix}
\title{
Accessor Method for the \code{PWM} slot in a
\code{\link{genomicProfiles}} object
}
\description{
Accessor Method for the \code{PWM} slot in a
\code{\link{genomicProfiles}} object
}
\usage{
PositionWeightMatrix(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfiles}}
    }
}
\details{
After creating a \code{\link{genomicProfiles}} object,
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
Patrick C.N. Martin and Nicolae Radu Zabe (2020) Dissecting the binding mechanisms of transcription factors to DNA using a statistical thermodynamics framework. CSBJ, 18, 3590-3605.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>
}

\examples{
#Loading data
data(ChIPanalyserData)
#Loading PFM files
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BEAF-32.pfm")
#Building genomicProfiles object
GPP<-genomicProfiles(PFM=PFM,PFMFormat="raw")
# Accessing Slot
PositionWeightMatrix(GPP)
}
