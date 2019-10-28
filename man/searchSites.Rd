\name{searchSites}
\alias{searchSites}


\title{
Searching function for Sites above threshold and predicted ChIP-seq Profiles
}
\description{
\code{searchSites} is function enabling quick extraction and search for
parameter combinations and/or loci in any \code{\link{genomicProfiles}} object
from computeOccupancy onwards.
}
\usage{
searchSites(Sites,lambdaPWM="all",BoundMolecules="all", Locus="all")
}

\arguments{
    \item{Sites}{
    \code{Sites} is either a \code{\link{genomicProfiles}} or the result of \code{computeOptimal}
    }
    \item{lambdaPWM}{
    \code{lambdaPWM} is a numeric vector describing the ScalingFactors
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
When testing numerous combinations of lambdaPWM and boundMolecules on top of many loci, it can
become challenging to navigate the large data output
\code{searchSites} will make searching in this slot a lot easier.
If all arguments are left at their default value of "all", then all Parameters
will be searched thus returning the full list of Sites above
threshold. If a value for lambdaPWM is user provided then only this lambdaPWM will be selected (all boundMolecules and loci will also be selected).
\code{searchSites} also works on the result of \code{computeOptimal}.
}
\value{
Returns object of same time as parsed to this function with only the parameters and/or loci selected.
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
  # Building genomicProfiles object
  GPP <- genomicProfiles(PFM=PFM, BPFrequency=DNASequenceSet)


  # Computing Genome Wide
  GenomeWide <- computeGenomeWideScore(genomicProfiles = GPP,
      DNASequenceSet = DNASequenceSet)

  #Compute PWM Scores
  PWMScores <- computePWMScore(genomicProfiles = GenomeWide,
       DNASequenceSet = DNASequenceSet, loci = eveLocus, chromatinState = Access)
  #Compute Occupnacy
  Occupancy <- computeOccupancy(genomicProfiles = PWMScores)
  searchSites(Occupancy,ScalingFactor=c(1,4), BoundMolecules = c(1,100),
      Locus="eve")

  #Compute ChIP profiles
  chipProfile <- computeChIPProfile(genomicProfiles=Occupancy,loci=eveLocus)
  searchSites(chipProfile,ScalingFactor=c(1,4), BoundMolecules = c(1,100),
      Locus="eve")

optimalParam <- computeOptimal(genomicProfiles = GPP,
        DNASequenceSet = DNASequenceSet,
        ChIPScore = eveLocusChip,
        chromatinState = Access,
        parameterOptions = OPP,
        parameter = "all",
        peakMethod="moving_kernel")

searchSites(optimalParam,ScalingFactor=c(1,4), BoundMolecules = c(1,100),
    Locus="eve")

}
}
