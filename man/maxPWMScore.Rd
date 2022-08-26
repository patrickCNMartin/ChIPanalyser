\name{maxPWMScore}


\alias{maxPWMScore}


\title{
Accessor function for \code{maxPWMScore} slot in a
\code{\link{genomicProfiles}} object.
}
\description{
Accessor function for \code{maxPWMScore} slot in a
\code{\link{genomicProfiles}} object.
}
\usage{
maxPWMScore(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfiles}} object.
    }
}
\details{
\code{maxPWMScore} is a numerical value that can be described as the
highest PWM score computed at a genome wide scale.
This value is computed and updated in the
\code{\link{genomicProfiles}} object after using the
\code{\link{computeGenomeWideScores}}.
}
\value{
Returns the value of assigned to the \code{maxPWMScore} slot in a
\code{\link{genomicProfiles}} object.
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
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BEAF-32.pfm")
#As an example of genome, this example will run on the Drosophila genome
\donttest{
if(!require("BSgenome.Dmelanogaster.UCSC.dm6", character.only = TRUE)){
    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6")
}
library(BSgenome.Dmelanogaster.UCSC.dm6)
DNASequenceSet <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6)
#Building data objects
GPP <- genomicProfiles(PFM=PFM,PFMFormat="JASPAR")

# Computing Genome Wide
GenomeWide <- computeGenomeWideScores(DNASequenceSet = DNASequenceSet,
    genomicProfiles = GPP)
maxPWMScore(GenomeWide)
## If used before computeGenomeWidePWMScore, will return NULL
}
}
