\name{ChIPanalyser-package}

\alias{ChIPanalyser}
\docType{package}
\title{
\packageTitle{ChIPanalyser}
}
\description{
\packageDescription{ChIPanalyser}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{ChIPanalyser}
\packageIndices{ChIPanalyser}

}
\author{
Patrick C.N. Martin <pm16057@essex.ac.uk>

And

Nicolae Radu Zabet <nzabet@essex.ac.uk>

Maintainer: \packageMaintainer{ChIPanalyser}
}


\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.

}

\keyword{ package }

\examples{
\donttest{
#Data extraction
data(ChIPanalyserData)
# path to Position Frequency Matrix
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BCDSlx.pfm")
#As an example of genome, this example will run on the Drosophila genome

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
    occupancy = Occupancy, occupancyProfileParameters = OPP)
#Estimating accuracy estimate
AccuracyEstimate <- profileAccuracyEstimate(LocusProfile = eveLocusChip,
    predictedProfile = chipProfile, occupancyProfileParameters = OPP)

    }
}
