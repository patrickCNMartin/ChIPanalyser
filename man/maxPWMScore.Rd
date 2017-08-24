\name{maxPWMScore}
\alias{maxPWMScore}
\alias{maxPWMScore-methods}


\title{
Accessor function for \code{maxPWMScore} slot in a
\code{\link{genomicProfileParameters}} object.
}
\description{
Accessor function for \code{maxPWMScore} slot in a
\code{\link{genomicProfileParameters}} object.
}
\usage{
maxPWMScore(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfileParameters}} object.
    }
}
\details{
\code{maxPWMScore} is a numerical value that can be described as the
highest PWM score computed at a genome wide scale.
This value is computed and updated in the
\code{\link{genomicProfileParameters}} object after using the
\code{\link{computeGenomeWidePWMScore}}.
}
\value{
Returns the value of assigned to the \code{maxPWMScore} slot in a
\code{\link{genomicProfileParameters}} object.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.

}
\author{
Patrick C. N. Martin <p.martin@essex.ac.uk>
}
\examples{

# Loading data
#Data extraction
data(ChIPanalyserData)
# path to Position Frequency Matrix
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BCDSlx.pfm")
#As an example of genome, this example will run on the Drosophila genome
\donttest{
if(!require("BSgenome.Dmelanogaster.UCSC.dm3", character.only = TRUE)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Dmelanogaster.UCSC.dm3")
}
library(BSgenome.Dmelanogaster.UCSC.dm3)
DNASequenceSet <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3)
#Building data objects
GPP <- genomicProfileParameters(PFM=PFM)

# Computing Genome Wide
GenomeWide <- computeGenomeWidePWMScore(DNASequenceSet = DNASequenceSet,
    genomicProfileParameters = GPP)
maxPWMScore(GenomeWide)
## If used before computeGenomeWidePWMScore, will return NULL
}
}
