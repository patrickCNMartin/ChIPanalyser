\name{averageExpPWMScore}
\alias{averageExpPWMScore}
\alias{averageExpPWMScore-methods}

\title{
Accessor for \code{averageExpPWMScore} slot in a
\code{\link{genomicProfileParameters}} object.
}
\description{
Extract or Access \code{averageExpPWMScore} slot in a
\code{\link{genomicProfileParameters}}
}
\usage{
averageExpPWMScore(object)
}
\arguments{
    \item{object}{
    \code{object} is a \code{genomicProfileParameters}
    }
}
\details{
As a general rule, \code{averageExpPWMScore} is computed and updated
internally by \code{\link{computeGenomeWidePWMScore}}.
Idealy, this slot should not be updated by user.
The \code{averageExpPWMScore} is the sum of the exponential of every PWM score
for a given DNA sequence and divided by the length of the said DNA sequence
(\code{\link{DNASequenceLength}}). This can either be the full
length sequence or only the accessible sequence
(see \code{\link{computeGenomeWidePWMScore}}).

}
\value{
Returns the \code{averageExpPWMScore} of a
\code{\link{genomicProfileParameters}} when computed.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.}
\author{
Patrick C.N. Martin <pm16057@essex.ac.uk>
}

\examples{
# Accessing Data
data(ChIPanalyserData)
# path to Position Frequency Matrix
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BCDSlx.pfm")
# Building genomicProfileParameters object
GPP <- genomicProfileParameters(PFM=PFM)
# Extracting AllSitesAboveThreshold slot
averageExpPWMScore(GPP)

## Note this slot is now empty as nothing has yet been computed
}
