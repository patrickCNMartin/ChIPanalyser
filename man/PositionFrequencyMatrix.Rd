\name{PositionFrequencyMatrix}
\alias{PositionFrequencyMatrix}
\alias{PositionFrequencyMatrix-methods}

\title{
Accessor method for the \code{PFM} slot in a
\code{genomicProfileParameters} object
}
\description{
Accessor method for the \code{PFM} slot in a
\code{genomicProfileParameters} object
}
\usage{
PositionFrequencyMatrix(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfileParameters}} object
    }
}
\details{
After creating a \code{\link{genomicProfileParameters}} object,
it is possible to access the Position Frequency Matrix slot.
However this slot will be empty if the \code{\link{genomicProfileParameters}}
object was built using directly a Position Weight Matrix.
See \code{\link{genomicProfileParameters}}
}
\value{
Returns the Position Frequency Matrix (\code{PFM} slot) used to compute the
\code{\link{PositionWeightMatrix}} in a
\code{\link{genomicProfileParameters}} object
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
PositionFrequencyMatrix(GPP)
}
