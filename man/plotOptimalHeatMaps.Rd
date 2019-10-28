\name{plotOptimalHeatMaps}
\alias{plotOptimalHeatMaps}

\title{
Heat Map of optimal Parameters
}
\description{
\code{plotOptimalHeatMaps} will plot heat maps of optimal
Parameters and highlight the optimal combination of
\code{\link{lambdaPWM}} and \code{\link{boundMolecules}}
}

\usage{
plotOptimalHeatMaps(optimalParam,contour=TRUE,col=NULL,main=NULL,layout=TRUE,overlay=FALSE)
}
\arguments{
    \item{optimalParam}{
    \code{optimalParam} is a list containing containing optimal matricies
    (or only one if only one paramter was selected). These matricies are the
    result of the \code{computeOptimal} function
    }
    \item{contour}{
    \code{parameter} is logical. Should contour lines be plotted?
    }
    \item{col}{
    \code{col} vector of colours to be used for each heat map.
    If none are specified, rainbow colours will be used. NOTE: colour
    vector will be recyled if not enough colours are provided.
    }
    \item{main}{\code{main} title.}
    \item{layout}{\code{layout} is either TRUE or FALSE specifying if stadard
    layout should be used or not. If TRUE, each heat map will be plotted on
    an individual page with a heat map scale of the right side.}
    \item{overlay}{\code{overlay} is either TRUE or FALSE specifying if an
    overlay plot should be produced. The overlay plot takes the top 10 percent
    of best performing parameters per scoring metric and overlays them in a
    single plot. The resulting plots shows the optimal set of paramters for all
    metrics combined.}
}
\details{
Once the optimal set of Parameters ( \code{\link{lambdaPWM}}
and \code{\link{boundMolecules}} ), it is possible to plot the results
in the form of a heat map. Each heat map will be plotted in a seperate page if
layout = TRUE, If layout= FALSE, it is up to the user to define how they wish
to layout there heat maps.
}
\value{
Returns a heat map of optimal combinations of \code{\link{lambdaPWM}}
and \code{\link{boundMolecules}}. The x axis represents the different
value assigned to lambda ( \code{\link{lambdaPWM}} )
and the y axis represents the different values to boundMolecules
( \code{\link{boundMolecules}} ).
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.}
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

#Building data objects
GPP <- genomicProfiles(PFM=PFM,BPFrequency=DNASequenceSet)

#Computing Optimal set of Parameters
optimalParam <- computeOptimal(genomicProfiles = GPP,
    DNASequenceSet = DNASequenceSet,
    ChIPScore = eveLocusChip,
    chromatinState = Access,
    parameterOptions = OPP,
    parameter = "all",
    peakMethod="moving_kernel")
plotOptimalHeatMaps(optimalParam)
}
}
