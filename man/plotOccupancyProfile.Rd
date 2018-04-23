\name{plotOccupancyProfile}
\alias{plotOccupancyProfile}

\title{
Plot Occupancy Profiles
}
\description{
\code{plotOccupancyProfile} plots the predicted profiles.
If provided, this functions will also plot ChIP-seq profiles,
PWMScores, DNAAccessibility, accuracy estimates and gene information.
}
\usage{
plotOccupancyProfile(predictedProfile,setSequence,
    chipProfile = NULL,DNAAccessibility = NULL, occupancy = NULL,PWM=FALSE,
    occupancyProfileParameters = NULL,geneRef = NULL,axis=TRUE,...)
}

\arguments{
    \item{predictedProfile}{
    \code{predictedProfile} is a \code{\link[GenomicRanges]{GRanges}}
    containing predicted ChIP-seq scores for ONE loci and ONE
    combination of \code{lambda}
    ( See \code{\link{ScalingFactorPWM}} ) and \code{\link{boundMolecules}}.
    }
    \item{setSequence}{
    \code{setSequence} is a \code{\link[GenomicRanges]{GRanges}} containing the
    locus of interest. This locus MUST be named.
    }

    \item{chipProfile}{
    \code{chipProfile} is a numeric vector containing real ChIP-seq scores.
    }
    \item{DNAAccessibility}{
    \code{DNAAccessibility} is a \code{\link[GenomicRanges]{GRanges}}
    containing accesible DNA sites
    }

    \item{occupancy}{
    \code{occupancy} is a \code{\link[GenomicRanges]{GRanges}} containing PWM
    Score and predicted Occupancies.
    }

    \item{PWM}{
    \code{PWM} is a logical value that in the case \code{occupancy} is
    provided which of occupancy scores of PWM scores hsould be plotted.
    Default set at FALSE
    }
    \item{occupancyProfileParameters}{
    \code{\link{occupancyProfileParameters}} is
    \code{\link{occupancyProfileParameters}} object containing parameters
    used for previous computations (see vignette for pipeline/function order)
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
it is possible to plot these profiles. The important aspect to keep in mind
with \code{plotOccupancyProfile} is that in only computes one profile
at a time. This means that the output of previous functions should be
subsetted carefully. This leaves to the user more flexibility in what
and how profiles should be plotted. Keep in mind that the output of
different function are more or less nested.
}
\value{
Returns a porfile plot with "Occupancy" on the y axis and DNA position on
the the X- axis. If the orange line is the predicted profile, the grey
shaded area represents real ChIP-seq data. The yellow boxes represent
regions on NON-accessible DNA. On the lower part of the plot, gene
information is plotted with respect to the strand they are localised on.
Finally, the blue vertical lines represent sites of either high Occupancy
or high PWM score depending on which had been selected. The minimal plot
will only contain the predicted profile.
The more data is provided the more will be plotted.
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
    occupancy = Occupancy, occupancyProfileParameters = OPP)


#Plotting Profile
plotOccupancyProfile(predictedProfile=chipProfile[[1]][[1]],
    setSequence=eveLocus,
    chipProfile = eveLocusChip[[1]],
    DNAAccessibility = Access,
    occupancy = AllSitesAboveThreshold(Occupancy)[[1]][[1]],
    occupancyProfileParameters = OPP,
    geneRef =geneRef)

## Changing graphical Parameters

#### In this examples ####
### predictedProfile will be reduce
### chipProfile will be blue
### and DNAAccessibility will be green
plotOccupancyProfile(predictedProfile=chipProfile[[1]][[1]],
    setSequence=eveLocus,
    chipProfile = eveLocusChip[[1]],
    DNAAccessibility = Access,
    occupancy = AllSitesAboveThreshold(Occupancy)[[1]][[1]],
    occupancyProfileParameters = OPP,
    geneRef =geneRef,col=c("red","blue","green"))

### If name is specified
plotOccupancyProfile(predictedProfile=chipProfile[[1]][[1]],
    setSequence=eveLocus,
    chipProfile = eveLocusChip[[1]],
    DNAAccessibility = Access,
    occupancy = AllSitesAboveThreshold(Occupancy)[[1]][[1]],
    occupancyProfileParameters = OPP,
    geneRef =geneRef,col=c("DNAAccessibility"="red","blue","green"))

}
}
