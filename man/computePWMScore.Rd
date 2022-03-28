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
computePWMScore(genomicProfiles,DNASequenceSet,
    loci = NULL, chromatinState = NULL,parameterOptions=NULL,cores=1, verbose = TRUE)
}

\arguments{
    \item{DNASequenceSet}{
    \code{DNASequenceSet} is a \code{\link[Biostrings]{DNAStringSet}} or a
    \code{BSgenome} containing the full sequence of the organism of interest.
    }
    \item{genomicProfiles}{
    \code{genomicProfiles} is a \code{\link{genomicProfiles}}
    object resulting from the \code{\link{computeGenomeWideScores}} function.
    }
    \item{loci}{
    \code{loci} is a \code{\link[GenomicRanges]{GRanges}}
    object containing the Loci of interest or a \code{\linkS4class{ChIPScore}} object
    result of \code{\link{processingChIP}} function.
    }
    \item{parameterOptions}{
    \code{parameterOptions} is a \code{\link{parameterOptions}}
    object containing parameters that you wish to parse/change when computing
    PWMScores.
    }
    \item{chromatinState}{
    \code{chromatinState} is a \code{\link[GenomicRanges]{GRanges}}
    object sites of accessible DNA or DNA affinity scores.
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
\code{\link{genomicProfiles}} object.
The \code{\link{profiles}} slot will have been updated.
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
if(!require("BSgenome.Dmelanogaster.UCSC.dm6", character.only = TRUE)){
    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6")
    }
library(BSgenome.Dmelanogaster.UCSC.dm6)
DNASequenceSet <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6)
eveLocusChip<-processingChIP(eveLocusChip,eveLocus)
#Building data objects
GPP <- genomicProfiles(PFM=PFM,BPFrequency=DNASequenceSet)


# Computing Genome Wide
GenomeWide <- computeGenomeWideScores(DNASequenceSet = DNASequenceSet,
    genomicProfiles = GPP)

#Compute PWM Scores
PWMScores <- computePWMScore(DNASequenceSet = DNASequenceSet,
    genomicProfiles = GenomeWide,
    loci  = eveLocusChip, chromatinState = Access)
PWMScores
}
}
