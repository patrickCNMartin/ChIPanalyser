\name{DNASequenceLength<-}
\alias{DNASequenceLength<-}
\alias{DNASequenceLength<--methods}

\title{
Setter Method for \code{DNASequenceLength} slot in a
\code{\link{genomicProfileParameters}} object.
}
\description{
Setter Method for \code{DNASequenceLength} slot in a
\code{\link{genomicProfileParameters}} object.}
\usage{
DNASequenceLength(object)<- value
}
\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfileParameters}} object.
    }
    \item{value}{
    \code{value} is a numerical value indicating the length of the DNA
    used and that will be assigned to \code{DNASequenceLength} slot in a
    \code{\link{genomicProfileParameters}}.
    }
}
\details{
Although it is possible to manually assign the value of the DNA Length,
it is not necessary as this value is computed internally in
\code{\link{computeGenomeWidePWMScore}}.
}
\value{
Returns a \code{\link{genomicProfileParameters}} object with an
updated value for the \code{DNASequenceLength} slot.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}
\author{
Patrick C.N. Martin <pm16057@essex.ac.uk>
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

# Building genomicProfileParameters object
GPP <- genomicProfileParameters(PFM=PFM,BPFrequency=DNASequenceSet)
# Setting new value for slot
DNASequenceLength(GPP)<- 1500000
# !! This value is computed internally !!
}
}
