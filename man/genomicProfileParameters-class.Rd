\name{genomicProfileParameters-class}
\Rdversion{1.1}
\docType{class}
\alias{genomicProfileParameters-class}




\alias{AllSitesAboveThreshold,genomicProfileParameters-method}

\alias{averageExpPWMScore,genomicProfileParameters-method}
\alias{BPFrequency<-,genomicProfileParameters,vector-method}
\alias{BPFrequency<-,genomicProfileParameters,DNAStringSet-method}

\alias{BPFrequency,genomicProfileParameters-method}
\alias{DNASequenceLength<-,genomicProfileParameters,vector-method}
\alias{DNASequenceLength,genomicProfileParameters-method}

\alias{maxPWMScore,genomicProfileParameters-method}

\alias{minPWMScore,genomicProfileParameters-method}
\alias{naturalLog<-,genomicProfileParameters,logical-method}
\alias{naturalLog,genomicProfileParameters-method}

\alias{NoAccess,genomicProfileParameters-method}
\alias{noOfSites<-,genomicProfileParameters,vector-method}
\alias{noOfSites,genomicProfileParameters-method}
\alias{PositionFrequencyMatrix<-,genomicProfileParameters,matrix-method}
\alias{PositionFrequencyMatrix<-,genomicProfileParameters,character-method}
\alias{PositionFrequencyMatrix,genomicProfileParameters-method}
\alias{PositionWeightMatrix<-,genomicProfileParameters,matrix-method}
\alias{PositionWeightMatrix,genomicProfileParameters-method}
\alias{PWMpseudocount<-,genomicProfileParameters,numeric-method}
\alias{PWMpseudocount,genomicProfileParameters-method}
\alias{PWMThreshold<-,genomicProfileParameters,numeric-method}
\alias{PWMThreshold,genomicProfileParameters-method}
\alias{ScalingFactorPWM<-,genomicProfileParameters,vector-method}
\alias{ScalingFactorPWM,genomicProfileParameters-method}

\alias{strandRule<-,genomicProfileParameters,character-method}
\alias{strandRule,genomicProfileParameters-method}
\alias{whichstrand<-,genomicProfileParameters,character-method}
\alias{whichstrand,genomicProfileParameters-method}
\alias{.ZeroBackground<-,genomicProfileParameters,vector-method}
\alias{.ZeroBackground,genomicProfileParameters-method}

\title{Class \code{"genomicProfileParameters"}}
\description{
\code{genomicProfileParameters} is a data/parameter storing object.
The main purpose of this object is to store parameters neccessary for
further computation but also store results. This makes parsing data and
parameters easier between the different of functions available in
\code{ChIPanalyser}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form
\code{genomicProfileParameters(PWM, PFM, ScalingFactorPWM, PWMpseudocount,
    BPFrequency, naturalLog, noOfSites, minPWMScore,
    maxPWMScore, PWMThreshold, AllSitesAboveThreshold,
    DNASequenceLength , averageExpPWMScore , strandRule ,
    whichstrand, NoAccess)}. \code{PWM} and \code{PFM} slots are the only
    slots that are absolutely required (At least one of the two).
    All other slots have a provided default value.
    We strongly advise you to tailor these parameters to increase goodness
    of fit between the predcited model and real ChIP-seq data.

}
\section{Slots}{
    \describe{
    \item{\code{PWM}:}{Object of class \code{"matrix"}:
    A Position Weight Matrix (either supplied or internally computed if
    PFM is provided)}
    \item{\code{PFM}:}{Object of class \code{"matrix"}:
    A Position Frequency Matrix (may also be a path to file containing PFM) }
    \item{\code{ScalingFactorPWM}:}{Object of class \code{"vector"}:
    A vector (or single value) contaning values for the ScalingFactorPWM
    (Also known as lambda).Default:1}
    \item{\code{PWMpseudocount}:}{Object of class \code{"numeric"}:
    A numeric value describing a PWMpseudocount for PWM computation.
    Default:1 }
    \item{\code{BPFrequency}:}{Object of class \code{"vector"}:
    Base Pair Frequency in the genome (if a DNA sequence is provided
    (as a \code{\link[Biostrings]{DNAStringSet}} or \code{BSgenome}),
    will be automatically computed internally). Default:c(0.25,0.25,0.25,0.25)}
    \item{\code{naturalLog}:}{Object of class \code{"logical"}:
    A logical value describing if natural Log will be used to compute
    the PWM (if FALSE then log2 will be used). Default: TRUE }
    \item{\code{noOfSites}:}{Object of class \code{"vector"}:
    A Positive integer descibing number of sites (in base pair) should
    be used from the PFM to compute PWM. Default =0 (Full width of
    binding site will be used when set to 0) }
    \item{\code{minPWMScore}:}{Object of class \code{"vector"}:
    Lowest PWM score accros the genome (computed and updated internally) }
    \item{\code{maxPWMScore}:}{Object of class \code{"vector"}:
    Highest PWM score across the genome (computed and updated internally) }
    \item{\code{PWMThreshold}:}{Object of class \code{"numeric"}:
    Threshold at which PWM Score should be selected (only sites above
    threshold will be selected - between 0 and 1) }
    \item{\code{AllSitesAboveThreshold}:}{Object of class \code{"GRList"}:
    Containins GRanges with sites above threshold and associated metrics
    (PWMscore and Occupancy) - Computed Internally}
    \item{\code{DNASequenceLength}:}{Object of class \code{"vector"}:
    Length of the Genome (or accesible genome) - computed internally }
    \item{\code{averageExpPWMScore}:}{Object of class \code{"vector"}:
    Average exponential PWM score across the genome
    (or accesible genome) - computed internally }
    \item{\code{strandRule}:}{Object of class \code{"character"}:
    "mean", "max" or "sum"  will dertermine how strand should be handle
    for computing PWM Scores. Default : "max" }
    \item{\code{whichstrand}:}{Object of class \code{"character"}:
    "+","-" or "+-" on which strand should PWM Score be computed.
    Default: "+-" }
    \item{\code{NoAccess}:}{Object of class \code{"vector"}:
    Stores Loci that do contain accesible DNA if it were to be the case
    (computed and updated internally) }
    \item{\code{ZeroBackground}:}{Object of class \code{"vector"}:
    Internal background value (computed internally)}

    }
}
\section{Extends}{
Class \code{\link{occupancyProfileParameters}}, directly.
}
\section{Methods}{
    \describe{

    \item{AllSitesAboveThreshold}{\code{signature(object =
    "genomicProfileParameters")}: Slot Accessor Method }

    \item{averageExpPWMScore}{\code{signature(object =
    "genomicProfileParameters")}: Slot Accessor Method}
    \item{BPFrequency<-}{\code{signature(object =
    "genomicProfileParameters", value = "vector")}:Slot Setter Method}
    \item{BPFrequency<-}{\code{signature(object =
    "genomicProfileParameters", value = "DNAStringSet")}:Slot Setter Method}

    \item{BPFrequency}{\code{signature(object =
    "genomicProfileParameters")}: Slot Accessor Method }
    \item{DNASequenceLength<-}{\code{signature(object =
    "genomicProfileParameters", value = "vector")}: Slot Setter Method}
    \item{DNASequenceLength}{\code{signature(object =
    "genomicProfileParameters")}: Slot Accessor Method }

    \item{maxPWMScore}{\code{signature(object =
    "genomicProfileParameters")}: Slot Accessor Method}

    \item{minPWMScore}{\code{signature(object =
    "genomicProfileParameters")}: Slot Accessor Method}
    \item{naturalLog<-}{\code{signature(object =
    "genomicProfileParameters", value = "logical")}: Slot Setter Method}
    \item{naturalLog}{\code{signature(object =
    "genomicProfileParameters")}: Slot Accessor Method }

    \item{NoAccess}{\code{signature(object =
    "genomicProfileParameters")}: Slot Accessor Method }
    \item{noOfSites<-}{\code{signature(object =
    "genomicProfileParameters", value = "vector")}: Slot Setter Method }
    \item{noOfSites}{\code{signature(object =
    "genomicProfileParameters")}: Slot Accessor Method}
    \item{PositionFrequencyMatrix<-}{\code{signature(object =
    "genomicProfileParameters", value = "matrix")}: Slot Setter Method }
    \item{PositionFrequencyMatrix<-}{\code{signature(object =
    "genomicProfileParameters", value = "character")}: Slot Setter Method }
    \item{PositionFrequencyMatrix}{\code{signature(object =
    "genomicProfileParam"}: Slot Accessor Method}
    \item{PositionWeightMatrix<-}{\code{signature(object =
    "genomicProfileParameters", value = "matrix")}:Slot Setter Method }
    \item{PositionWeightMatrix}{\code{signature(object =
    "genomicProfileParameters")}: Slot Accessor Method}
    \item{PWMpseudocount<-}{\code{signature(object =
    "genomicProfileParameters", value = "numeric")}: Slot Setter Method}
    \item{PWMpseudocount}{\code{signature(object =
    "genomicProfileParameters")}:Slot Accessor Method }
    \item{PWMThreshold<-}{\code{signature(object =
    "genomicProfileParameters", value = "numeric")}: Slot Setter Method }
    \item{PWMThreshold}{\code{signature(object =
    "genomicProfileParameters")}: Slot Accessor Method}
    \item{ScalingFactorPWM<-}{\code{signature(object =
    "genomicProfileParameters", value = "vector")}: Slot Setter Method }
    \item{ScalingFactorPWM}{\code{signature(object =
    "genomicProfileParameters")}: Slot Accessor Method }
    \item{strandRule<-}{\code{signature(object =
    "genomicProfileParameters", value = "character")}: Slot Setter Method }
    \item{strandRule}{\code{signature(object =
    "genomicProfileParameters")}: Slot Accessor Method }
    \item{whichstrand<-}{\code{signature(object =
    "genomicProfileParameters", value = "character")}: Slot Setter Method }
    \item{whichstrand}{\code{signature(object =
    "genomicProfileParameters")}: Slot Accessor Method }

    }
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.

}
\author{
Partick C. N. Martin <p.martin@essex.ac.uk>
}

\seealso{
\code{\link{genomicProfileParameters}}
\code{\link{occupancyProfileParameters-class}}
\code{\link{occupancyProfileParameters}}


}
\examples{

#Data extraction
data(ChIPanalyserData)
# path to Position Frequency Matrix
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BCDSlx.pfm")
#As an example of genome, this example will run on the Drosophila genome

if(!require("BSgenome.Dmelanogaster.UCSC.dm3", character.only = TRUE)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Dmelanogaster.UCSC.dm3")
}
library(BSgenome.Dmelanogaster.UCSC.dm3)
DNASequenceSet <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3)

# Building genomicProfileParameters object
GPP <- genomicProfileParameters(PFM=PFM, ScalingFactorPWM=c(1,1.5,2),
    BPFrequency=DNASequenceSet)
}
\keyword{classes}
