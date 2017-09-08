\name{genomicProfileParameters}
\alias{genomicProfileParameters}


\title{
Genomic Profile Parameter Object
}
\description{
\code{genomicProfileParameter} is a data/parameter storing object.
The main purpose of this object is to store parameters neccessary for further
computation but also store results. This makes parsing data and parameters
easier between the different of functions available in \code{ChIPanalyser}.

}
\usage{


genomicProfileParameters(PWM = NULL, PFM = NULL,PFMFormat="raw",
    ScalingFactorPWM = 1,
    PWMpseudocount = 1, noOfSites = 0, BPFrequency = rep(0.25, 4),
    naturalLog = FALSE, PWMThreshold = 0.7, strandRule = "max",
    whichstrand = "+-")

}

\arguments{
    \item{PFM}{ \code{PFM} is a positon Frequency matrix.
    This matrix can either be supplied directly as a matrix
    (in the following form: 1 row for each base pair - ACTG in order -
    and a number of rows equal to the number of base pairs in the binding site)
    or can be supplied as path to a file containing the
    Position Frequency Matrix. At least ONE of either the
    PFM or PWM must be supplied. This is the minimum to build a
    \code{genomicProfileParameters} object.
    NOTE: Although the object will build, NO computation can be carried
    out without at least a PFM or a PWM.

    }
    \item{PWM}{ \code{PWM} is a position weight matrix.
    This can be provided in matrix form (1 row for each base pair -
    ACTG in order - and a number of rows equal to the number of base pairs
    in the binding site).
    If a PFM is provided the PWM will be computed internally.
    NOTE: Although the object will build, NO computation can be carried
    out without at least a PFM or a PWM.

    }
    \item{PFMFormat}{
    \code{PFMFormat} is a character string describing the format of the
    file used (if any) in PositionFrequencyMatrix. ChIPanalyser
    handles most PFM formats thus this argument should take one of the
    following: raw, transfac, JASPAR or sequences.
    }
    \item{ScalingFactorPWM}{ \code{\link{ScalingFactorPWM}}:
    A vector (or single value) contaning values for the ScalingFactorPWM
    (Also known as lambda).Default:1
    }
    \item{PWMpseudocount}{ \code{\link{PWMpseudocount}}:
    A numeric value describing a PWMpseudocount for PWM computation. Default:1
    }
    \item{BPFrequency}{ \code{\link{BPFrequency}}: Base Pair Frequency in the
    genome (if a DNA sequence is provided
    (as a \code{\link[Biostrings]{DNAStringSet}} or \code{BSgenome}) ,
    will be automatically computed internally). Default:c(0.25,0.25,0.25,0.25
    }
    \item{naturalLog}{ \code{\link{naturalLog}}:
    A logical value describing if natural Log will be used to compute
    the PWM (if FALSE then log2 will be used). Default: TRUE
    }
    \item{noOfSites}{ \code{\link{noOfSites}}:
    A Positive integer descibing number of sites (in base pair) should
    be used from the PFM to compute PWM. Default =0 (Full width of binding
    site will be used when set to 0)
    }

    \item{PWMThreshold}{ \code{\link{PWMThreshold}}:
    Threshold at which PWM Score should be selected
    (only sites above threshold will be selected - between 0 and 1)
    }

    \item{strandRule}{ \code{\link{strandRule}}:
    `mean', `max' or `sum'  will dertermine how strand should be
    handle for computing PWM Scores. Default : `max'
    }
    \item{whichstrand}{ \code{\link{whichstrand}}: `+',`-' or `+-'
    on which strand should PWM Score be computed. Default: `+-'
    }

}
\details{
Most functions in \code{ChIPanalyser} require some form of a
\code{genomicProfileParameters} object. It is required to build a
genomicProfileParameters object before starting any computation.
Most functions will return a genomicProfileParameters and thus this
object should be parsed to the next function and so on. Most parameters
are customisable and we strongly advise to do so. This can either be
done directly (advised method) when building the
\code{genomicProfileParameters} or by using setter methods.
The former is advised as certain value are computed internally and
require that parameters are set before hand ( e.g \code{BPFrequency} ).
It is important to note that not all slots/parameters can be set by user.
The value assinged to these slots are computed internally.
For more information on this slots, see
\code{\link{genomicProfileParameters-class}}.
The following describes the default values and which slot are
computed internally:

PWM = To be supplied by user


PFM = To be supplied by user


PFMFormat ="raw"


ScalingFactorPWM = 1


PWMpseudocount = 1


BPFrequency = c(0.25,0.25,0.25,0.25)


naturalLog = TRUE


noOfSites = 0


minPWMScore = Internally Computed


maxPWMScore = Internally Computed


PWMThreshold = 0.7


AllSitesAboveThreshold = Internally Computed


DNASequenceLength = Internally Computed (may still be provided by user see
\code{\link{DNASequenceLength<-}})


averageExpPWMScore = Internally Computed


strandRule = `max'


whichstrand = `+-'


NoAccess = Internally Computed


ZeroBackground = Internally Computed



}
\value{
Returns a \code{genomicProfileParameters} object
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.

}
\author{
Patrick C. N. Martin  <p.martin@essex.ac.uk>
}





\seealso{

\code{\link{genomicProfileParameters-class}}
\code{\link{occupancyProfileParameters}}
\code{\link{occupancyProfileParameters-class}}
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
