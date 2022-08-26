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
    occupancy = NULL,goodnessOfFit = NULL,PWM=FALSE,
    geneRef = NULL,addLegend = TRUE,...)
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
    \item{addLegend}{
    \code{addLegend} is a logical value defining if the legend should be added.
    The legend will add all elements provided. See details. 
    }
    
    \item{...}{
    Any other graphical Parameter of the following : cex, cex.lab, cex.main, 
    densityCS , densityGR , ylab, xlab, main, colPred,
    colChIP, colOccup, colCS, colGR, n_axis_ticks. See details. 
    }

}
\details{
Once the predicted ChIP-seq like profiles have been computed,
it is possible to plot these profiles.

This functions allows to control graphical parameters. In short:

* col = color values - exact number of colors or colors that will 
be used in a colorRampPalettte.

* cex = font sizes - for text, axis labels and main 

* Density = fill density for chromatin state and/or geneRef blocks

Pred = predictedProfile
ChIP = ChIP score (Experimental ChIP data)
CS = Chromatin States
GR = Gene reference 
Occup = Occupnacy locations
}
\value{
Returns a profile plot with "Occupancy" on the y axis and DNA position on
the the X- axis.

}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
Patrick C.N. Martin and Nicolae Radu Zabe (2020) 
Dissecting the binding mechanisms of transcription factors to DNA using a statistical thermodynamics framework. CSBJ, 18, 3590-3605.
}
\author{
Patrick C.N. Martin <pcnmartin@gmail.com>
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
#Building data objects
GPP <- genomicProfiles(PFM=PFM,PFMFormat="JASPAR", BPFrequency=DNASequenceSet)


# Computing Genome Wide
GenomeWide <- computeGenomeWideScores(DNASequenceSet = DNASequenceSet,
    genomicProfiles = GPP)

#Compute PWM Scores
PWMScores <- computePWMScore(DNASequenceSet = DNASequenceSet,
    genomicProfiles = GenomeWide,
    loci = top, chromatinState = Access)
#Compute Occupnacy
Occupancy <- computeOccupancy(genomicProfiles = PWMScores)

#Compute ChIP profiles
chipProfile <- computeChipProfile(loci = top,
    genomicProfiles = Occupancy)


#Plotting Profile
plotOccupancyProfile(predictedProfile=chipProfile,
    ChIPScore = chip,
    chromatinState = Access,
    occupancy = Occupancy,
    geneRef =geneRef)

plotOccupancyProfile(predictedProfile=chipProfile,
    ChIPScore = chip,
    chromatinState = Access,
    occupancy = Occupancy,
    geneRef = geneRef,
    colCS = c("red","blue"),
    densityGR = 60)


}
}
