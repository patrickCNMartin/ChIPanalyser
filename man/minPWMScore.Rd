\name{minPWMScore}
\alias{minPWMScore}


\title{
Accessor method the \code{minPWMScore} slot in a
\code{\link{genomicProfiles}} object
}
\description{
Accessor method the \code{minPWMScore} slot in a
\code{\link{genomicProfiles}} object
}
\usage{
minPWMScore(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfiles}} object.
    }
}
\details{
\code{minPWMScore} can be described as the lowest PWM score computed at
a genome wide scale. Although it is possible to assigne a value
to \code{minPWMScore}, we strongly advise to use the value
computed and assigned internally. This value is computed in the
\code{\link{computeGenomeWideScores}} function.
}
\value{
Returns the value assigned to the \code{minPWMScore} slot in a
\code{\link{genomicProfiles}} object.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.}
\author{
Patrick C. N. Martin <p.martin@essex.ac.uk>
}
\examples{

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
GPP <- genomicProfiles(PFM=PFM,PFMFormat="raw")

# Computing Genome Wide
GenomceWide <- computeGenomeWideScores(DNASequenceSet = DNASequenceSet,
    genomicProfiles = GPP)
minPWMScore(GenomceWide)

## If used before computeGenomeWidePWMScore, will return NULL

}
}
