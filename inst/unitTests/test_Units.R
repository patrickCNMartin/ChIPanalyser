## Test Units ##

## Loading test data##

library("ChIPanalyser")

data(ChIPanalyserData)
# path to Position Frequency Matrix
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BEAF-32.pfm")
#As an example of genome, this example will run on the Drosophila genome

if(!require("BSgenome.Dmelanogaster.UCSC.dm3", character.only = TRUE)){
    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm3")
}
library(BSgenome.Dmelanogaster.UCSC.dm3)
DNASequenceSet <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3)

load(file.path(path.package("ChIPanalyser"), "unitTests", "GPPWithPFM.rda"))
# S4 Class testing

#test_genomicProfileParametersBuild <- function(){
  #  GPP <- genomicProfileParameters()
  #  checkTrue(class(GPP)=="genomicProfileParameters",
#    #"Class Build genomicProfileParameters: OK")
#}

#test_occupancyProfileParametersBuild <- function(){
#    OPP <- occupancyProfileParameters()
#    checkTrue(class(OPP)=="occupancyProfileParameters",
#    "Class Build occupancyProfileParameters: OK")
#}

# PFM to PWM in genomicProfileParameters object

#test_GPPInternal <- function(){
#    load(file.path(path.package("ChIPanalyser"), "unitTests",
#    "GPPWithPFM.rda"))
#    GPP <- genomicProfileParameters(PFM=PFM,BPFrequency=DNASequenceSet)
#    checkIdentical(PositionFrequencyMatrix(GPPWithPFM),
#    PositionFrequencyMatrix(GPP))
#    checkIdentical(PositionWeightMatrix(GPPWithPFM),
#    PositionWeightMatrix(GPP))
#    checkEquals(BPFrequency(GPPWithPFM),BPFrequency(GPP))
#}

# Testing Genome Wide scoring No Access

#test_GenomeWideScoreNoAccess <- function(){
  #  load(file.path(path.package("ChIPanalyser"), "unitTests",
  #  "genomeWideScoringNoAccess.rda"))
  #  GPP <- genomicProfileParameters(PFM=PFM,BPFrequency=DNASequenceSet)
  #  GW <- computeGenomeWidePWMScore(DNASequenceSet,GPP)
  #  checkTrue(class(GW)=="genomicProfileParameters",
  #  "Class Check genomicProfileParameters: OK")
  #  #checkEquals(maxPWMScore(genomeWideScoringNoAccess),maxPWMScore(GW))
  #  #checkEquals(minPWMScore(genomeWideScoringNoAccess),minPWMScore(GW))
    #checkEquals(averageExpPWMScore(genomeWideScoringNoAccess),
    #averageExpPWMScore(GW))

#}


# Testing Genome Wide Scoring with Access

#test_GenomeWideScoreAccess <- function(){
#    load(file.path(path.package("ChIPanalyser"), "unitTests",
#    "genomeWideScoringAccess.rda"))
#    GPP <- genomicProfileParameters(PFM=PFM,BPFrequency=DNASequenceSet)
#    GW <- computeGenomeWidePWMScore(DNASequenceSet,GPP,Access,GenomeWide=FALSE)
#    checkTrue(class(GW)=="genomicProfileParameters",
#    "Class Check genomicProfileParameters: OK")
    #checkEquals(maxPWMScore(genomeWideScoringAccess),maxPWMScore(GW))
    #checkEquals(minPWMScore(genomeWideScoringAcess),minPWMScore(GW))
    #checkEquals(averageExpPWMScore(genomeWideScoringAccess),
    #averageExpPWMScore(GW))

#}


# Testing Sites above Threshold No Access

#test_SitesAboveThresholdNoAccess <- function(){
#    load(file.path(path.package("ChIPanalyser"), "unitTests",
#    "PWMScoreNoAccess.rda"))
#    GPP <- genomicProfileParameters(PFM=PFM,BPFrequency=DNASequenceSet)
#    GW <- computeGenomeWidePWMScore(DNASequenceSet,GPP)
#    SAT <- computePWMScore(DNASequenceSet,GPP,top)
#    checkTrue(class(SAT)=="genomicProfileParameters")

    #checkTrue(is(AllSitesAboveThreshold(SAT), "GRangesList"),
    #"Class Check genomicProfileParameters: OK")
    #checkEquals(AllSitesAboveThreshold(SAT)$PWMScores,
    #AllSitesAboveThreshold(PWMScoresNoAccess)[[1]]$PWMScores)
    #checkEquals(start(AllSitesAboveThreshold(SAT)[[1]]),
    #AllSitesAboveThreshold(PWMScoresNoAccess)[[1]]$PWMScores)
#}

# Testing sites above Threshold with Access
#test_SitesAboveThresholdAccess <- function(){
#    load(file.path(path.package("ChIPanalyser"), "unitTests",
#    "PWMScoreAccess.rda"))
#    GPP <- genomicProfileParameters(PFM=PFM,BPFrequency=DNASequenceSet)
#    GW <- computeGenomeWidePWMScore(DNASequenceSet,GPP,Access,GenomeWide=FALSE)
#    SAT <- computePWMScore(DNASequenceSet,GPP,top,Access)
#    checkTrue(class(SAT)=="genomicProfileParameters")

    #checkTrue(is(AllSitesAboveThreshold(SAT), "GRangesList"),
    #"Class Check genomicProfileParameters: OK")
    #checkEquals(AllSitesAboveThreshold(SAT)[[1]]$PWMScores,
    #AllSitesAboveThreshold(PWMScoresAccess)[[1]]$PWMScores)
    #checkEquals(start(AllSitesAboveThreshold(SAT)[[1]]),
    #start(AllSitesAboveThreshold(PWMScoresAccess)[[1]]))
#}

# Testing Occupancy with No access

#test_OccupancyNoAccess <- function(){
#    load(file.path(path.package("ChIPanalyser"), "unitTests",
#    "OccupancyNoAccess.rda"))
#    GPP <- genomicProfileParameters(PFM=PFM,BPFrequency=DNASequenceSet)
#    GW <- computeGenomeWidePWMScore(DNASequenceSet,GPP)
#    SAT <- computePWMScore(DNASequenceSet,GPP,top)
#    OPP <- occupancyProfileParameters()
#    Occup <- computeOccupancy(SAT,OPP)
#    checkTrue(class(Occup)=="genomicProfileParameters")

    #checkEquals(AllSitesAboveThreshold(Occup)[[1]]$Occupancy,
    #AllSitesAboveThreshold(OccupancyNoAccess)[[1]]$Occupancy)

#}

# Testing Occupancy with Access
#test_OccupancyAccess <- function(){
#    load(file.path(path.package("ChIPanalyser"), "unitTests",
#    "OccupancyAccess.rda"))
#    GPP <- genomicProfileParameters(PFM=PFM,BPFrequency=DNASequenceSet)
#    GW <- computeGenomeWidePWMScore(DNASequenceSet,GPP)
#    SAT <- computePWMScore(DNASequenceSet,GPP,top)
#    OPP <- occupancyProfileParameters()
#    Occup <- computeOccupancy(SAT,OPP)
#    checkTrue(class(Occup)=="genomicProfileParameters")

    #checkEquals(AllSitesAboveThreshold(Occup)[[1]]$Occupancy,
    #AllSitesAboveThreshold(OccupancyAccess)[[1]]$Occupancy)

#}
