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
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BEAF-32.pfm")
#As an example of genome, this example will run on the Drosophila genome

if(!require("BSgenome.Dmelanogaster.UCSC.dm6", character.only = TRUE)){
    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6")
}
library(BSgenome.Dmelanogaster.UCSC.dm6)
DNASequenceSet <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6)
#Building data objects
GPP <- genomicProfiles(PFM=PFM,PFMFormat="raw",BPFrequency=DNASequenceSet)

chip<-processingChIP(chip,top)
# Computing Genome Wide
GenomeWide <- computeGenomeWideScores(DNASequenceSet = DNASequenceSet,
    genomicsProfiles = GPP)

#Compute PWM Scores
PWMScores <- computePWMScore(genomicsProfiles = GenomeWide,
    DNASequenceSet = DNASequenceSet,
    loci = top, chromatinState = Access)
#Compute Occupnacy
Occupancy <- computeOccupancy(genomicsProfiles = PWMScores,
    parameterOptions = OPP)

#Compute ChIP profiles
chipProfile <- computeChIPProfile(genomicProfiles = Occupancy,
     loci = top,
    parameterOptions = OPP)
#Estimating accuracy estimate
AccuracyEstimate <- profileAccuracyEstimate(genomicProfiles = chipProfile,
     ChIPScore = chip,
     parameterOptions = OPP)

    }
}
