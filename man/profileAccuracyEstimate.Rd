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
    occupancyProfileParameters = NULL,method="all")
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
    \item{method}{\code{method} is the method that will be used to assess model
    quality agianst ChIP-seq data. Method can be one of the following:
    pearson, spearman, kendall, ks, geometric,fscore or all. All of these will
    also return the MSE. Fscore contains f-score, precision,recall, MCC,
    Accuracy and AUC ROC.}
}
\details{
In order to assess the quality of the model against experimental ChIP-seq data,
ChIPanalyser offers a wide range of method to choose from. These methods are
also used when computing optimal paramters.
}
\value{
Returns a list of two elements. The first element represents lists containing
the model quality assessments for every combination of parameters (Bound Molecules and lambda)
for every genomic region. The second element of the list contains the result of the ROCR package:
False positives, False Negative, etc... See ROCR package.
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
    source("https://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Dmelanogaster.UCSC.dm3")
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
