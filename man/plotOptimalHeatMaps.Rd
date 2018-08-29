\name{plotOptimalHeatMaps}
\alias{plotOptimalHeatMaps}

\title{
Heat Map of optimal Parameters
}
\description{
\code{plotOptimalHeatMaps} will plot heat maps of optimal
Parameters and highlight the optimal combination of
\code{\link{ScalingFactorPWM}} and \code{\link{boundMolecules}}
}

\usage{
plotOptimalHeatMaps(optimalParam, parameter = "all", Contour = TRUE)
}
\arguments{
    \item{optimalParam}{
    \code{optimalParam} is a list resulting of \code{\link{computeOptimal}}.
    This list should contain three elements: optimal Parameters,
    optimal Matrix and parameter.
    }
    \item{parameter}{
    \code{parameter} is a character string representing which parameter
    should be plotted. There are four possibilities:
    "correlation","MSE","theta" and "all"
    }
    \item{Contour}{
    \code{Contour} is a logical value indicating if contour
    lines should be plotted.
    }
}
\details{
Once the optimal set of Parameters ( \code{\link{ScalingFactorPWM}}
and \code{\link{boundMolecules}} ), it is possible to plot the results
in the form of a heat map. There are four possible heat maps:
a "correlation" heat map that will only show the correlation between a
predicted Profile and true ChIP-seq profiles for each combination of
Parameters ; "MSE" will show a similar heat map only with the
Mean Squared Error associated to each predicted profile,
true profile and parameter combination; "theta" is a in house metric
describing the best fit between high correlation and low MSE
(see \code{\link{thetaThreshold}}); finally "all" will plot all of the above.
}
\value{
Returns a heat map of optimal combinations of \code{\link{ScalingFactorPWM}}
and \code{\link{boundMolecules}}. The x axis represents the different
value assigned to lambda ( \code{\link{ScalingFactorPWM}} )
and the y axis represents the different values to boundMolecules
( \code{\link{boundMolecules}} ). The hilighted box describes the
optimal combination of parameters to either maximise "correlation"
and "theta" or to minimise "MSE".
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
GPP <- genomicProfileParameters(PFM=PFM,BPFrequency=DNASequenceSet)
OPP <- occupancyProfileParameters()
#Computing Optimal set of Parameters
optimalParam <- computeOptimal(DNASequenceSet = DNASequenceSet,
    genomicProfileParameters = GPP,
    LocusProfile = eveLocusChip,
    setSequence = eveLocus,
    DNAAccessibility = Access,
    occupancyProfileParameters = OPP,
    parameter = "all")
plotOptimalHeatMaps(optimalParam, parameter="all")
}
}
