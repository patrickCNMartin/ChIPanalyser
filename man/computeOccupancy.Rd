\name{computeOccupancy}
\alias{computeOccupancy}

\title{
Compute Occupancy values from PWM Scores based on model.
}
\description{
\code{computeOccupancy} will compute the Occupancy from PWM Scores.
As described in detail in the vignette,
\code{ChIPanalyser} uses PWM Scores, DNA Accessibility data, the number of
bound molecules and a sclaing factor of Transcription Factor specificty.
This function will compute occupancy using the
values assigned to each variable.
}
\usage{
computeOccupancy(genomicProfiles,parameterOptions = NULL,
    norm = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{genomicProfiles}{
    \code{genomicProfiles} is a a \code{\link{genomicProfiles}}
    object resulting from \code{\link{computePWMScore}}.
    IT is important to use this resulting object as the occuapancy will only be
    computed for sites above a threshold.
    }
    \item{parameterOptions}{
    \code{parameterOptions} is a
    \code{\link{parameterOptions}} object containing the adequate
    values assigned to each Parameter. If not Supplied
    (\code{parameterOptions} = NULL), a new object will be created
    internally using default values.
    }
    \item{norm}{

    \code{norm} a logical value which determines
    if the occupancy should be normalised or not.
    }
    \item{verbose}{
    \code{verbose} a logical value which determines
    if progress messages are printed or not.
    }
}
\details{
\code{computeOccupancy} will compute the Occupancy from PWM Scores.
As described in detail in the vignette,
\code{ChIPanalyser} uses PWM Scores, DNA Accessibility data,
the number of bound molecules and a sclaing factor of
Transcription Factor specificty.
This function will compute occupancy using the values assigned
to each variable. It should also be noted that the
\code{\link{parameterOptions}}
object contains a set of parameters used to compute Occupancy
(not only restricted to this ). These parameters are often dependant on
real ChIP-Seq data and will influence
the goodness of fit between the predicted model an real ChIP-seq data.
We strongly advise that the values assigned to each parameter should be
customiszed in order to increase the model ageement with
real world biological data.

}
\value{
\code{computeOccupancy} will return a \code{\link{genomicProfiles}}.
The main difference will reside in the
\code{\link{profiles}} slot.
This slot is generally a list or \code{\link[GenomicRanges]{GRangesList}}.
Within these list type structures are enclosed
\code{\link[GenomicRanges]{GRanges}} containing the positions of site
above threshold, PWMScores and Occupancy for each site.
The series of GRanges  will depend on the number of loci that are
tested and the number of element in the list will depend on the various
combinations of \code{lambdaPWM} and \code{boundMolecules}.


}
\references{

Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
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
GPP <- genomicProfiles(PFM=PFM,PFMFormat="JASPAR",BPFrequency=DNASequenceSet)
OPP <- parameterOptions()
# Computing Genome Wide
GenomeWide <- computeGenomeWideScores(genomicProfiles = GPP,
    DNASequenceSet = DNASequenceSet)

#Compute PWM Scores
PWMScores <- computePWMScore(genomicProfiles = GenomeWide,
    DNASequenceSet = DNASequenceSet,
    loci = top,
    chromatinState = Access)
#Compute Occupnacy
Occupancy <- computeOccupancy(genomicProfiles = PWMScores,
    parameterOptions = OPP)
Occupancy

}
}
