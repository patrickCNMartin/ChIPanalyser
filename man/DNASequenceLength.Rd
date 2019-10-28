\name{DNASequenceLength}


\alias{DNASequenceLength}

\title{
Accessor method for \code{DNASequenceLength} slot in a
\code{\link{genomicProfiles}}
}
\description{
Accessor method for \code{DNASequenceLength} slot in a
\code{\link{genomicProfiles}}

}
\usage{
DNASequenceLength(object)
}

\arguments{
    \item{object}{
    \code{object} is  a \code{\link{genomicProfiles}}
    }
}
\details{
The model on which is based \code{ChIPanalyser} requires the length of the
DNA sequence used to compute scores. In this circustance,
this DNA Length is the total length of the DNA of the organism of interest
or the the Accessible DNA at a genome wide scale.
}
\value{
Returns \code{DNASequenceLength} slot in a
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

#Data extraction
data(ChIPanalyserData)
# path to Position Frequency Matrix
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BCDSlx.pfm")
#As an example of genome, this example will run on the Drosophila genome
\donttest{
if(!require("BSgenome.Dmelanogaster.UCSC.dm3", character.only = TRUE)){
    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm3")
}
library(BSgenome.Dmelanogaster.UCSC.dm3)
DNASequenceSet <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3)

#Building genomicProfiles object
GPP <- genomicProfiles(PFM=PFM,BPFrequency=DNASequenceSet)
# Computing Genome Wide
GenomceWide <- computeGenomeWideScores(DNASequenceSet = DNASequenceSet,
    genomicProfiles = GPP)

DNASequenceLength(GenomceWide)

}
}
