\name{plotOccupancyProfile}
\alias{plotOccupancyProfile}

\title{
Plot Occupancy Profiles
}
\description{
\code{plotOccupancyProfile} plots the predicted profiles.
If provided, this functions will also plot ChIP-seq profiles,
PWMScores (or Occupancy), chromatin States, Goodness of Fit estimates and gene information.
}
\usage{
  plotOccupancyProfile(predictedProfile, ChIPScore = NULL,chromatinState = NULL,
      asCS = FALSE,occupancy = NULL,goodnessOfFit = NULL,PWM=FALSE,
      geneRef = NULL,axis=TRUE,...)
}

\arguments{
    \item{predictedProfile}{
    \code{predictedProfile} is a either \code{\link[GenomicRanges]{GRanges}}
     containing the predicted profiles for one loci, all loci selected for one paramter,
     or all loci selected for all parameter combinations selected. (see \code{\link{searchSites}})
    }
    \item{ChIPScore}{
    \code{ChIPScore} is a \code{ChIPscore} object containing ChIPscore (or
    a list of numeric values representing ChIP scores (Experimental ChIP))
    }
    \item{chromatinState}{
    \code{chromatinState} is a \code{\link[GenomicRanges]{GRanges}}
    containing accesible DNA sites or chromatin States.
    }
    \item{asCS}{
    \code{asCS} is a logical defining if chromatinStates should be plotted as
    accessibility data (Default FALSE) or if chromatin states should be plotted
    in seperate track.
    }

    \item{occupancy}{
    \code{occupancy} is a \code{\link[GenomicRanges]{GRanges}} or a \code{genomicProfiles}
    object contaning PWM scores and Occupancy ( see \code{computeOccupancy})
    }
    \item{goodnessOfFit}{ \code{goodnessOfFit} results of the
    \code{\link{profileAccuracyEstimate}} function.}

    \item{PWM}{
    \code{PWM} is a logical value that in the case \code{occupancy} is
    provided which of occupancy scores of PWM scores hsould be plotted.
    Default set at FALSE
    }

    \item{geneRef}{
    \code{geneRef} is a \code{\link[GenomicRanges]{GRanges}}
    containing gene information on exons,introns, UTR's, enhancers or any
    other genetic element to be plotted.
    }
    \item{axis}{
    \code{axis} is a logical value indicating if axes should be plotted or not.
    }
    \item{...}{
    Any other graphical Parameter of the following : col, density, border,
    lty, lwd, cex, cex.axis, xlab, ylab, xlim, ylim, las and axislables.
    Each Parameter will be parsed in the same order as the arguments to
    this function. If the name of the argument is specified, the argument value
    will be parsed to the correct internal plotting. See examples below
    }

}
\details{
Once the predicted ChIP-seq like profiles have been computed,
it is possible to plot these profiles.
}
\value{
Returns a profile plot with "Occupancy" on the y axis and DNA position on
the the X- axis. If the orange line is the predicted profile, the grey
shaded area represents real ChIP-seq data. The yellow boxes represent
regions on NON-accessible DNA. On the lower part of the plot, gene
information is plotted with respect to the strand they are localised on.
Finally, the blue vertical lines represent sites of either high Occupancy
or high PWM score depending on which had been selected. The minimal plot
will only contain the predicted profile.
The more data is provided the more will be plotted.
IMPORTANT: the colours described above represent the default colour values.
If colours are provided by the user, they will be updated in plotting order or
if the argument has been speicified assigned to said argument.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.}
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
if(!require("BSgenome.Dmelanogaster.UCSC.dm6", character.only = TRUE)){
    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6")
    }
library(BSgenome.Dmelanogaster.UCSC.dm6)
DNASequenceSet <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6)
#Building data objects
GPP <- genomicProfiles(PFM=PFM,BPFrequency=DNASequenceSet)


# Computing Genome Wide
GenomeWide <- computeGenomeWideScores(DNASequenceSet = DNASequenceSet,
    genomicProfiles = GPP)

#Compute PWM Scores
PWMScores <- computePWMScore(DNASequenceSet = DNASequenceSet,
    genomicProfiles = GenomeWide,
    loci = eveLocus, chromatinState = Access)
#Compute Occupnacy
Occupancy <- computeOccupancy(genomicProfiles = PWMScores)

#Compute ChIP profiles
chipProfile <- computeChipProfile(loci = eveLocus,
    genomicProfiles = Occupancy)


#Plotting Profile
plotOccupancyProfile(predictedProfile=chipProfile,
    ChIPScore = eveLocusChip,
    chromatinState = Access,
    occupancy = Occupancy,
    geneRef =geneRef)

## Changing graphical Parameters

#### In this examples ####
### predictedProfile will be red
### chipProfile will be blue
### and DNAAccessibility will be green
plotOccupancyProfile(predictedProfile=chipProfile,
    ChIPScore = eveLocusChip,
    chromatinState = Access,
    occupancy = Occupancy,
    geneRef =geneRef,col=c("red","blue","green"))

### If name is specified
plotOccupancyProfile(predictedProfile=chipProfile,
    ChIPScore = eveLocusChip,
    chromatinState = Access,
    occupancy = Occupancy,
    geneRef =geneRef,col=c("DNAAccessibility"="red","blue","green"))

}
}
