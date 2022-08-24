\name{computeGenomeWideScores}
\alias{computeGenomeWideScores}

\title{
Computing Genome Wide scores
}
\description{
\code{computeGenomeWideScores} compute the max and min PWM score over the entire genome.
}
\usage{
computeGenomeWideScores(genomicProfiles, DNASequenceSet, chromatinState = NULL, parameterOptions = NULL, cores = 1, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genomicProfiles}{
  \code{genomicProfiles} is a \code{genomicProfiles} object containing the PFM, PWM of interest.
}
  \item{DNASequenceSet}{
  \code{DNASequenceSet} is a BSgenome or DNAStringSet containing the sequence of the organism of interest.
}
  \item{chromatinState}{
  \code{chromatinState} is a GRanges object containing the chromatin States. This can either represent regions of accessible DNA or Chromatin state affinities.
}
  \item{parameterOptions}{
   \code{parameterOptions} is a \code{parameterOptions} object containing parameters that you wish to change. The genomicProfiles object will be updated using the values assigned to \code{parameterOptions}
}
  \item{cores}{
 \code{cores} is the number or cores that will be used (Numeric value - Default = 1 )
}
  \item{verbose}{
  \code{verbose} is a logical value that will determine if internal progress message will be printed.
}
}
\details{
\code{computeGenomeWideScores} function computes PWM scores over the entire genome (or accessible Genome if chromatin State are provided ). Genome wide scores are used to determine the maximum and minimum PWM score as well as the average exponential score. These scores will in turn be used to determine which score are above the PWM theshold. The average exponential score is an integrale part of the equation used to compute Occupancy. Using defualt settings, ChIPanalyser will only compute occupancy on the top 70\% of PWM scores. This threshold can be changed. See \code{\link{PWMThreshold}}
}
\value{
Returns a \code{genomicsProfiles} object with updated values for max score, min score and averageExpPWMScore.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
Patrick C.N. Martin and Nicolae Radu Zabe (2020) Dissecting the binding mechanisms of transcription factors to DNA using a statistical thermodynamics framework. CSBJ, 18, 3590-3605.
}
\author{
Patrick C.N Martin <pm16057@essex.ac.uk>
}

\examples{

  \donttest{
  if(!require("BSgenome.Dmelanogaster.UCSC.dm6", character.only = TRUE)){
      if (!requireNamespace("BiocManager", quietly=TRUE))
          install.packages("BiocManager")
      BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6")
      }
  library(BSgenome.Dmelanogaster.UCSC.dm6)
  DNASequenceSet <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6)
  # Building genomicProfiles object
  GPP <- genomicProfiles(PFM=PFM, BPFrequency=DNASequenceSet)


  # Computing Genome Wide
  GenomeWide <- computeGenomeWideScores(genomicProfiles = GPP,
      DNASequenceSet = DNASequenceSet)
  }

}
