\name{PositionFrequencyMatrix}


\alias{PositionFrequencyMatrix}

\title{
Accessor method for the \code{PFM} slot in a
\code{genomicProfiles} object
}
\description{
Accessor method for the \code{PFM} slot in a
\code{genomicProfiles} object
}
\usage{
PositionFrequencyMatrix(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfiles}} object
    }
}
\details{
After creating a \code{\link{genomicProfiles}} object,
it is possible to access the Position Frequency Matrix slot.
However this slot will be empty if the \code{\link{genomicProfiles}}
object was built using directly a Position Weight Matrix.
See \code{\link{genomicProfiles}}
}
\value{
Returns the Position Frequency Matrix (\code{PFM} slot) used to compute the
\code{\link{PositionWeightMatrix}} in a
\code{\link{genomicProfiles}} object
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
PositionFrequencyMatrix(GPP)
}
