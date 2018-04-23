\name{computePWMScore}
\alias{computePWMScore}

\title{
Compute PWM Scores of sites above threshold.
}
\description{
\code{computePWMScore} will compute and extract all sites that exhibit a
PWM Score higher than a threshold.
This threshold (see \code{\link{PWMThreshold}}) will determine the percentage
of total sites that should NOT be considered.
}
\usage{
computePWMScore(DNASequenceSet, genomicProfileParameters,
    setSequence = NULL, DNAAccessibility = NULL,cores=1, verbose = TRUE)
}

\arguments{
    \item{DNASequenceSet}{
    \code{DNASequenceSet} is a \code{\link[Biostrings]{DNAStringSet}} or a
    \code{BSgenome} containing the full sequence of the organism of interest.
    }
    \item{genomicProfileParameters}{
    \code{genomicProfileParameters} is a \code{\link{genomicProfileParameters}}
    object resulting from the \code{\link{computeGenomeWidePWMScore}} function.
    }
    \item{setSequence}{
    \code{setSequence} is a \code{\link[GenomicRanges]{GRanges}}
    object containing the Loci of interest.
    }
    \item{DNAAccessibility}{
    \code{DNAAccessibility} is a \code{\link[GenomicRanges]{GRanges}}
    object containing Accesible DNA.
    }
    \item{cores}{
    \code{cores} is the number of cores used to compute PWM Scores. 
    }
    \item{verbose}{
    \code{verbose} is a logical value indicating if progress messages
    should be printed or not.
    }
}
\details{
After determining genome wide scores, it is possible to only compute and
extract high affinity sites (in the sense that they have a high PWM Score).
If a \code{\link{PWMThreshold}} is not set by user,
the default value is set at 0.7.
This means that 70 \% of sites will NOT be selected.
Only the top 30 \% will be computed and extracted.
If one is interested in all PWM Scores at a genome wide scale
( or accessible DNA ), this is possible by setting
\code{\link{PWMThreshold}} to zero.
}
\value{
\code{\link{computePWMScore}} will return a
\code{\link{genomicProfileParameters}} object.
The \code{\link{AllSitesAboveThreshold}} slot will have been updated.
This slot will now contain a \code{\link[GenomicRanges]{GRangesList}}
with each element being a \code{\link[GenomicRanges]{GRanges}}.
This GRanges will contain postion of each sites
(start, end and strand) and the PWMScore associated to that site.
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
    source("https://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Dmelanogaster.UCSC.dm3")
    }
library(BSgenome.Dmelanogaster.UCSC.dm3)
DNASequenceSet <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3)

#Building data objects
GPP <- genomicProfileParameters(PFM=PFM,BPFrequency=DNASequenceSet)


# Computing Genome Wide
GenomeWide <- computeGenomeWidePWMScore(DNASequenceSet = DNASequenceSet,
    genomicProfileParameters = GPP)

#Compute PWM Scores
PWMScores <- computePWMScore(DNASequenceSet = DNASequenceSet,
    genomicProfileParameters = GenomeWide,
    setSequence = eveLocus, DNAAccessibility = Access)
PWMScores
}
}
