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
  profileAccuracyEstimate(genomicProfiles,ChIPScore,
      parameterOptions=NULL,method="all",cores=1)
}

\arguments{
    \item{genomicProfiles}{
    \code{genomicProfiles} is the result of \code{\link{computeChIPProfile}}
    }
    \item{ChIPScore}{
    \code{ChIPScore} is the result of \code{processingChIP}. Extracted/Normalised
    experimental ChIP scores.

    }
    \item{parameterOptions}{
    \code{parameterOptions} is a
    \code{\link{parameterOptions}} object for paramter specification.
    }
    \item{method}{\code{method} is the method that will be used to assess model
    quality agianst ChIP-seq data. Method can be one of the following:
    pearson, spearman, kendall, ks, geometric,fscore, MSE,or all.Fscore contains f-score, precision,recall, MCC,
    Accuracy and AUC ROC.}
    \item{cores}{
    \code{cores} is the number of cores used to extract ChIP scores. Default = 1
    }
}
\details{
In order to assess the quality of the model against experimental ChIP-seq data,
ChIPanalyser offers a wide range of method to choose from. These methods are
also used when computing optimal paramters.
}
\value{
Returns list of goodness of fit metrics for each loci and each parameter selected.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
Patrick C.N. Martin and Nicolae Radu Zabe (2020) Dissecting the binding mechanisms of transcription factors to DNA using a statistical thermodynamics framework. CSBJ, 18, 3590-3605.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>
}
\examples{

#Data extraction
data(ChIPanalyserData)
# path to Position Frequency Matrix
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BEAF-32.pfm")
#As an example of genome, this example will run on the Drosophila genome
\donttest{
  if(!require("BSgenome.Dmelanogaster.UCSC.dm6", character.only = TRUE)){
      if (!requireNamespace("BiocManager", quietly=TRUE))
          install.packages("BiocManager")
      BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6")
      }
  library(BSgenome.Dmelanogaster.UCSC.dm6)
  DNASequenceSet <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6)
  # Building genomicProfiles object
  GPP <- genomicProfiles(PFM=PFM,PFMFormat="JASPAR", BPFrequency=DNASequenceSet)


  # Computing Genome Wide
  GenomeWide <- computeGenomeWideScore(genomicProfiles = GPP,
      DNASequenceSet = DNASequenceSet)

  #Compute PWM Scores
  PWMScores <- computePWMScore(genomicProfiles = GenomeWide,
       DNASequenceSet = DNASequenceSet, loci = top, chromatinState = Access)
  #Compute Occupnacy
  Occupancy <- computeOccupancy(genomicProfiles = PWMScores)

  #Compute ChIP profiles
  chipProfile <- computeChIPProfile(genomicProfiles=Occupancy,loci=top)
#Estimating accuracy estimate
AccuracyEstimate <- profileAccuracyEstimate(genomicProfiles = chipProfile,
    ChIPScore = chip,
    occupancyProfileParameters = OPP)
}
}
