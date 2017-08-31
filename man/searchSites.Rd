\name{searchSites}
\alias{searchSites}


\title{
Searching function for Sites above threshold and predicted ChIP-seq Profiles
}
\description{
\code{searchSites} is function enabling quick extraction and search for
the \code{\link{genomicProfileParameters}}'s
\code{\link{AllSitesAboveThreshold}}
slot or  a list resulting from \code{\link{computeChipProfile}}.
}
\usage{
searchSites(Sites,ScalingFactor="all",BoundMolecules="all", Locus="all")
}

\arguments{
    \item{Sites}{
    \code{Sites} is either a \code{\link{genomicProfileParameters}} object
    or a list ( result of \code{\link{computeChipProfile}})
    }
    \item{ScalingFactor}{
    \code{ScalingFactor} is a numeric vector describing the ScalingFactors
    that should be searched within \code{Sites}.
    }
    \item{BoundMolecules}{
    \code{BoundMolecules} is a numeric vector describing the BoundMolecules
    that should be searched within \code{Sites}.
    }
    \item{Locus}{
    \code{Locus} is a character vector describing the Loci that should be
    searched within Sites.
    }
}
\details{
After PWM Scores and/or occupancy have been computed
(see \code{\link{computePWMScore}} and \code{\link{computeOccupancy}},
it is possible to extract and search the
\code{\link{AllSitesAboveThreshold}} slot. As this
slot may be fairly large it can become difficult to navigate.
\code{searchSites} will make searching in this slot a lot easier.
If all arguments are left at their default value of "all", then all Parameters
will be searched thus returning the full list of Sites above
threshold.
\code{searchSites} also works on the result of \code{computeChipProfile}.
}
\value{
Returns a list or \code{GRangesList}.
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
GPP <- genomicProfileParameters(PFM=PFM,BPFrequency=DNASequenceSet,
    ScalingFactorPWM=c(1,2,3,4))
OPP <- occupancyProfileParameters(boundMolecules=c(1,10,100))

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
# Search on occupancy
searchSites(Occupancy,ScalingFactor=c(1,4), BoundMolecules = c(1,100))
#Compute ChIP profiles
chipProfile <- computeChipProfile(setSequence = eveLocus,
    occupancy = Occupancy,
    occupancyProfileParameters = OPP)

## Search on chipProfile
searchSites(Occupancy,ScalingFactor=c(1,4), BoundMolecules = c(1,100),
    Locus="eve")

}
}
