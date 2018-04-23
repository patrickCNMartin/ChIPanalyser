\name{computeChipProfile}
\alias{computeChipProfile}

\title{
Computing ChIP-seq like profiles from Occupancy data.
}
\description{
\code{computeChipProfile} compute ChIP-seq like profile from occupancy data.
Occupancy data is computed using \code{\link{computeOccupancy}}.
}
\usage{
computeChipProfile(setSequence, occupancy, occupancyProfileParameters = NULL,
    norm = TRUE, method = c("moving_kernel","truncated_kernel","exact"),
    peakSignificantThreshold= NULL,cores=1, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{setSequence}{
    \code{setSequence} is a \code{\link[GenomicRanges]{GRanges}} containing the
    loci of interest.
    It is strongly advised to name each loci / range in the
    \code{\link[GenomicRanges]{GRanges}}.
    }
    \item{occupancy}{
    \code{occupancy} is a \code{link{genomicProfileParameters}} object result
    of the \code{\link{computeOccupancy}} function.
    This \code{\link{genomicProfileParameters}} object should contain a list,
    GRanges or a GRangesList with both PWMScore and Occupancy
    as meta data columns. To check if your object contains the right data,
    please see \code{\link{AllSitesAboveThreshold}}.
    }
    \item{occupancyProfileParameters}{
    \code{occupancyProfileParameters} is an
    \code{\link{occupancyProfileParameters}}
    object containing the desired values for each Parameter. If left as NULL,
    \code{computeChipProfile} will generate a new
    \code{occupancyProfileParameters}
    object using the default values
    (see \code{\link{occupancyProfileParameters}}).
    }
    \item{norm}{
    \code{norm} is a logical value. If TRUE, the ChIP-seq like profile will be
    normalised towards maximum Occupancy.
    If FALSE, the profile will be left as is.
    }
    \item{method}{
    \code{method} is a character string  of one of the following:
    c("moving_kernel","truncated_kernel","exact"). If set to moving_kernel,
    the peaks will be approximated using Rcpp (Default).
    If set to truncated_kernel, the peaks will be approximated
    however this method does not require Rcpp. If set to exact, the peaks will
    not be approximated.

    }
    \item{peakSignificantThreshold}{
    \code{peakSignificantThreshold} is a threshold at which peaks will be
    selected. IMPORTANT: if you select "moving_kernel" as described in
    \code{method} then this threshold is a numeric value describing the
    peak tail hight cutoff value (Default = 0.001).
    In the case of "truncated_kernel" and "exact", the threhsold represents
    a distance in base pair from the peak summit at which the peak
    should be cut (Default = 1250). The default is set to NULL in this
    function. This just means that either the value is provided bu user with
    the appropriate method. If not, the default will be selected depending
    on the method selected.
    }
    \item{cores}{
    \code{cores} is the number of cores that will be used to compute
    ChIP profiles.
    }
    \item{verbose}{
    \code{verbose} is a logical value. If TRUE, progress messages will be
    displayed in console.
    If FALSE, no progress messages will be dispalyed in console.
    }
}
\details{
\code{computeChipProfile} converts Transcription Factor occuapncy to a profile
resembling the one of a ChIP-seq profile. A certain set of Parameters are
required in order to build ChIP profiles.
These Parameters are defined and storedin a
\code{\link{occupancyProfileParameters}} object. These parameters are:
\code{\link{chipMean}}, \code{\link{chipSd}}, \code{\link{chipSmooth}},
\code{\link{stepSize}},\code{\link{backgroundSignal}},
\code{\link{maxSignal}} and \code{\link{removeBackground}}.
All these Parameters have default values already stored.
However, for an optimal fit, it is advised to derive these values
from actual ChIP-seq data.
For more information on these parameters,
see \code{\link{occupancyProfileParameters}}.
This functions also requires a set of sequencesin form of a
\code{\link[GenomicRanges]{GRanges}}. The sequence set are the loci of interest
on which the ChIP-seq profile will be computed.

}
\value{
Returns a list containing all ChIP-seq like profile for every combination of
\code{\link{ScalingFactorPWM}} and \code{\link{boundMolecules}}.
The correlation and Mean Squared Error between the prdicted ChIP profile
and actual ChIP-seq profile for the same loci
will vary depending on the value given for \code{\link{ScalingFactorPWM}}
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}
\author{
Patrick C.N. Martin <pm16057@essex.ac.uk>
}


\examples{

#Extracting Data
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
# Building genomicProfileParameters object
GPP <- genomicProfileParameters(PFM=PFM, BPFrequency=DNASequenceSet)

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
chipProfile
}
}
