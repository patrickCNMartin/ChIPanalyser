\name{profileAccuracyEstimate}
\alias{profileAccuracyEstimate}

\title{
Estimating Accuracy of predicted Profiles
}
\description{
\code{profileAccuracyEstimate} will compare the predicted ChIP-seq-like
profile to real ChIP-seq data and return a set of metrics describing how
accurate the predicted model is compared to real data.
}
\usage{
profileAccuracyEstimate(LocusProfile, predictedProfile,
    occupancyProfileParameters = NULL)
}

\arguments{
    \item{LocusProfile}{
    \code{LocusProfile} is a list containing normalised ChIP-seq profiles of
    the Loci of interest (real data). Each is profile is a numeric vector of
    length equals to the length of the locus in base pair.
    }
    \item{predictedProfile}{
    \code{predictedProfile} is the result of the
    \code{\link{computeChipProfile}} function.
    Generally, the output of this function comes in the form of a
    list of \code{\link[GenomicRanges]{GRangesList}}.
    Each GRangesList contains a \code{\link[GenomicRanges]{GRanges}}
    with the predicted ChIP-seq-like profiles for each Loci of interest.
    }
    \item{occupancyProfileParameters}{
    \code{occupancyProfileParameters} is a
    \code{\link{occupancyProfileParameters}} object
    }
}
\details{
The accuracy of the predicted profile may be estimated by measuring
corraltion, Mean Squared Error and theta (in house metric based on a
modified ratio of correlation over MSE) between predicted Profiles and
real ChIP-seq data. Actual ChIP-seq profiles should be normalised to a
base pair level (Enrichement divded by the width of the range for that
given score - the end result is a numeric vector of length equals to the
length of the locus in base pairs). It should be noted that if an
\code{\link{occupancyProfileParameters}} object is not supplied,
then one will be created internally. However, we strongly advise to use
the same occupancyProfileParameters object used previously.
}
\value{
Returns a list of lists. Each element in the list represents a combination
of lambda (see \code{\link{ScalingFactorPWM}}) and bound molecules
(see \code{\link{boundMolecules}}) and the list within each element is
he list of Loci of interest. Finally, at the core of these lists is a
named vector containing correlation and MSE for the given Loci but also
meanCorr, meanMSE and meanTheta for all loci for a given combination of
Lambda and bound molecules.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>
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
    genomicProfileParameters = GenomeWide,
    setSequence = eveLocus, DNAAccessibility = Access)
#Compute Occupnacy
Occupancy <- computeOccupancy(AllSitesPWMScore = PWMScores,
    occupancyProfileParameters = OPP)

#Compute ChIP profiles
chipProfile <- computeChipProfile(setSequence = eveLocus,
    occupancy = Occupancy,
    occupancyProfileParameters = OPP)
#Estimating accuracy estimate
AccuracyEstimate <- profileAccuracyEstimate(LocusProfile = eveLocusChip,
    predictedProfile = chipProfile,
    occupancyProfileParameters = OPP)
}
}
