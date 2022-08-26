\name{computeChIPProfile}
\alias{computeChIPProfile}

\title{
Computing ChIP-seq like profiles from Occupancy data.
}
\description{
\code{computeChIPProfile} compute ChIP-seq like profile from occupancy data.
Occupancy data is computed using \code{\link{computeOccupancy}}.
}
\usage{
computeChIPProfile(genomicProfiles, loci, parameterOptions = NULL,
    norm = TRUE, method = c("moving_kernel","truncated_kernel","exact"),
    peakSignificantThreshold= NULL,cores=1, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{genomicProfiles}{
    \code{genomicProfiles} is the result of \code{\link{computeOccupancy}}.
    This object should be a \code{\link{genomicProfiles}} object.
    }
    \item{loci}{
    \code{loci} is either a \code{\link{GRanges}} or \code{\linkS4class{ChIPScore}}
    object. ChIPScore-class will be the result of \code{\link{processingChIP}}.
    This object represents the set of Loci you are interested in analysing.
    If you have followed the full ChIPanalyser pipe line, you would have used
    the processingChIP function that would return a ChIPScore-class object containing
    your loci of interest. GRanges are also supported if you are only using
    part of the pipeline.
    }
    \item{parameterOptions}{
    \code{parameterOptions} is a \code{\link{parameterOptions}} object.
    This object is used to store the numerous paramters offered by ChIPanalyser.
    This argument is optional as all arguments are also parse in both ChIPScore-class
    and genomicProfiles objects. If you wanted to make some last minute changes,
    \code{\link{parameterOptions}} is the way to go. We recomend that you set
    your desired options before hand.
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
\code{computeChIPProfile} converts Transcription Factor occuapncy to a profile
resembling the one of a ChIP-seq profile. Internally a few paramters are required
to build a ChIP like profile. These parameters are either defined and stored in
a \code{\linkS4class{ChIPScore}} object (Paramters are updated based on
your ChIP data ), a \code{\link{genomicProfiles}} (user defined at the
start of the analysis) or a \code{\link{parameterOptions}} (if you want to
update values as you go along)
}
\value{
Returns a \code{\link{genomicProfiles}} objec containing all ChIP-seq like
profile for every combination of \code{\link{lambdaPWM}} and \code{\link{boundMolecules}} provided by the user.

}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
Patrick C.N. Martin and Nicolae Radu Zabe (2020) Dissecting the binding mechanisms of transcription factors to DNA using a statistical thermodynamics framework. CSBJ, 18, 3590-3605.
}
\author{
Patrick C.N. Martin <pcnmartin@gmail.com>
}


\examples{

#Extracting Data
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
GPP <- genomicProfiles(PFM=PFM, PFMFormat="JASPAR",BPFrequency=DNASequenceSet)


# Computing Genome Wide
GenomeWide <- computeGenomeWideScores(genomicProfiles = GPP,
    DNASequenceSet = DNASequenceSet)

#Compute PWM Scores
PWMScores <- computePWMScore(genomicProfiles = GenomeWide,
     DNASequenceSet = DNASequenceSet, loci = top, chromatinState = Access)
#Compute Occupnacy
Occupancy <- computeOccupancy(genomicProfiles = PWMScores)

#Compute ChIP profiles
chipProfile <- computeChIPProfile(genomicProfiles=Occupancy,loci=top)
chipProfile
}
}
