\name{AllSitesAboveThreshold}
\alias{AllSitesAboveThreshold}
\alias{AllSitesAboveThreshold-methods}


\title{
Accessor for \code{AllSitesAboveThreshold} slot in a
\code{\link{genomicProfileParameters}} object.
}
\description{
Extract or Access \code{AllSitesAboveThreshold} slot in
\code{\link{genomicProfileParameters}} object.
}
\usage{
AllSitesAboveThreshold(object)
}

\arguments{
    \item{object}{
    \code{\link{genomicProfileParameters}} object.
    }
}
\details{
As a general rule, \code{AllSitesAboveThreshold} is computed and updated
internally in \code{\link{computePWMScore}}, \code{\link{computeOccupancy}},
\code{\link{computeChipProfile}}. Idealy, this slot should not be updated
by user. There are three differenrent occurances of this slot.
As a \code{\link{GRanges}}, \code{AllSitesAboveThreshold} represents the sites
of \emph{one} locus that are above \code{\link{PWMThreshold}}.
As a \code{\link{GRangesList}}, \code{AllSitesAboveThreshold} represents the
sites of \emph{multiple} loci that are above \code{\link{PWMThreshold}}.
As a \code{list}, \code{AllSitesAboveThreshold} represents a list of
\code{\link{GRanges}} or \code{\link{GRangesList}} (as described above).
When computing occupancy \code{\link{computeOccupancy}}, it is possible to
use mutliple values of \code{\link{ScalingFactorPWM}} and
\code{\link{boundMolecules}}.
Each element in the list reprensent a combination of
\code{\link{ScalingFactorPWM}} and \code{\link{boundMolecules}}.
The specific combination of each element
can be found by using \code{\link{names}} of \code{AllSitesAboveThreshold}
or by using \code{\link{ScalingFactorPWM}} and \code{\link{boundMolecules}}.
}
\value{
Returns a \code{GRanges} , \code{GRangesList} or \code{list} containing
containing sites with a Position Weight Matrix Score
above a \code{\link{PWMThreshold}}.
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
AllSitesAboveThreshold(GPP)

## Note this slot is now empty as nothing has yet been computed
}
