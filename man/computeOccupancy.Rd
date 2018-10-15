\name{computeOccupancy}
\alias{computeOccupancy}

\title{
Compute Occupancy values from PWM Scores based on model.
}
\description{
\code{computeOccupancy} will compute the Occupancy from PWM Scores.
As described in detail in the vignette,
\code{ChIPanalyser} uses PWM Scores, DNA Accessibility data, the number of
bound molecules and a sclaing factor of Transcription Factor specificty.
This function will compute occupancy using the
values assigned to each variable.
}
\usage{
computeOccupancy(AllSitesPWMScore,DNAAccessibility=NULL,
    occupancyProfileParameters = NULL,
    norm = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{AllSitesPWMScore}{
    \code{AllSitesPWMScore} is a a \code{\link{genomicProfileParameters}}
    object resulting from \code{\link{computePWMScore}}.
    IT is important to use this resulting object as the occuapancy will only be
    computed for sites above a threshold.
    }
    \item{DNAAccessibility}{
    \code{DNAAccessibility} is a GRanges object containing Ranges of accesible
    sequences. This object ccould also contain a DNA affinity score as a
    metadata column.
    }
    \item{occupancyProfileParameters}{
    \code{occupancyProfileParameters} is a
    \code{\link{occupancyProfileParameters}} object containing the adequate
    values assigned to each Parameter. If not Supplied
    (\code{occupancyProfileParamters} = NULL),a new object will be created
    internally using default values.
    }
    \item{norm}{

    \code{norm} a logical value which determines
    if the occupancy should be normalised or not.
    }
    \item{verbose}{
    \code{verbose} a logical value which determines
    if progress messages are printed or not.
    }
}
\details{
\code{computeOccupancy} will compute the Occupancy from PWM Scores.
As described in detail in the vignette,
\code{ChIPanalyser} uses PWM Scores, DNA Accessibility data,
the number of bound molecules and a sclaing factor of
Transcription Factor specificty.
This function will compute occupancy using the values assigned
to each variable. It should also be noted that the
\code{\link{occupancyProfileParameters}}
object contains a set of parameters used to compute Occupancy
(not only restricted to this ). These parameters are often dependant on
real ChIP-Seq data and will influence
the goodness of fit between the predicted model an real ChIP-seq data.
We strongly advise that the values assigned to each parameter should be
customiszed in order to increase the model ageement with
real world biological data.

}
\value{
\code{computeOccupancy} will return a \code{\link{genomicProfileParameters}}.
The main difference will reside in the
\code{\link{AllSitesAboveThreshold}} slot.
This slot is generally a list or \code{\link[GenomicRanges]{GRangesList}}.
Within these list type structures are enclosed
\code{\link[GenomicRanges]{GRanges}} containing the positions of site
above threshold, PWMScores and Occupancy for each site.
The series of GRanges  will depend on the number of loci that are
tested and the number of element in the list will depend on the various
combinations of \code{ScalingFactorPWM} and \code{boundMolecules}.


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

#Building data objects
GPP <- genomicProfileParameters(PFM=PFM,BPFrequency=DNASequenceSet)
OPP <- occupancyProfileParameters()
# Computing Genome Wide
GenomeWide <- computeGenomeWidePWMScore(DNASequenceSet = DNASequenceSet,
    genomicProfileParameters = GPP)

#Compute PWM Scores
PWMScores <- computePWMScore(DNASequenceSet = DNASequenceSet,
    genomicProfileParameters = GenomeWide, setSequence = eveLocus,
    DNAAccessibility = Access)
#Compute Occupnacy
Occupancy <- computeOccupancy(AllSitesPWMScore = PWMScores,
    occupancyProfileParameters = OPP)
Occupancy

}
}
