\name{computeGenomeWidePWMScore}
\alias{computeGenomeWidePWMScore}

\title{
Computing genome wide PWM scores.
}
\description{
\code{computeGenomeWidePWMScore} scores the entire genome with a
Position Weight Matrix.
It is possible to restrict the scoring of the genome to only Accessible DNA.
}
\usage{
computeGenomeWidePWMScore(DNASequenceSet, genomicProfileParameters,
    DNAAccessibility = NULL, verbose = TRUE)
}

\arguments{
    \item{DNASequenceSet}{
    \code{DNASequenceSet} is a \code{BSgenome} or a
    \code{\link[Biostrings]{DNAStringSet}}.
    The DNASequenceSet is the sequence of the organism of interest.
    For a question of ease, it is advised to use a
    \code{\link[Biostrings]{DNAStringSet}}
    }
    \item{genomicProfileParameters}{
    \code{genomicProfileParameters} is a
    \code{\link{genomicProfileParameters}} object.
    This object contains parameters needed to compute
    \code{computeGenomeWidePWMScore}.
    GenomicProfileParameters also contain the Position Weight Matrix.
    See \code{\link{genomicProfileParameters}} for more information.
    }
    \item{DNAAccessibility}{
    \code{DNAAccessibility} is a \code{\link[GenomicRanges]{GRanges}} object
    containing Accessible DNA. If provided, the score will be score only on
    accesible DNA instead of genome wide.
    }

    \item{verbose}{
    \code{verbose} is a logical value. If TRUE, progress messages will be
    displayed in console. If FALSE, no progress message will be displayed.
    }
}
\details{
After creating a \code{\link{genomicProfileParameters}} object,
the genome can be scored in order to extract information
required for further analysis. These values will be \code{\link{maxPWMScore}},
\code{\link{minPWMScore}}, \code{\link{averageExpPWMScore}}
and finally  \code{\link{DNASequenceLength}}.
These values will be used in \code{\link{computePWMScore}}.
}
\value{
Returns a \code{\link{genomicProfileParameters}} object.
This object will contain updated value for certain slots
\code{\link{maxPWMScore}}, \code{\link{minPWMScore}},
\code{\link{averageExpPWMScore}} and finally
\code{\link{DNASequenceLength}}.
It should be mentioned that certain slot will remain empty as they will
be computed in other functions.
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
DNASequenceSet <-DNAStringSet(BSgenome.Dmelanogaster.UCSC.dm3[["chr2R"]])
names(DNASequenceSet) <- "chr2R"

# Building genomicProfileParameters object
GPP <- genomicProfileParameters(PFM=PFM, BPFrequency=DNASequenceSet)

# Computing Genome Wide
GenomeWide <- computeGenomeWidePWMScore(DNASequenceSet = DNASequenceSet,
    genomicProfileParameters = GPP)
# Computing Genome wide with DNA accessibility
GenomeWide <- computeGenomeWidePWMScore(DNASequenceSet = DNASequenceSet,
    genomicProfileParameters = GPP, DNAAccessibility = Access,
    GenomeWide = FALSE)

}
}
