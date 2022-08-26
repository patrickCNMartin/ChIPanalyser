\name{BPFrequency<-}

\alias{BPFrequency<-}




\title{
Setter method for \code{BPFrequency} slot in a
\code{\link{genomicProfiles}} object.
}
\description{
Setter method for \code{BPFrequency} slot in a
\code{\link{genomicProfiles}} object.
If base pair frequency is unknown, \code{BPFrequency} will compute base pair
frequency from a DNA sequence.
}

\usage{
BPFrequency(object)<-value
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfiles}} object.
    }
    \item{value}{
    \code{value} can three different objects:


    A vector of length 4 containing the probability of occurrence of each base
    pair (A,C,G,T) in order.
    Default value is c(0.25,0.25,0.25,0.25).


    A \code{BSgenome} of the organism of interest. The base pair frequency will
    automatically be computed and updated in
    \code{\link{genomicProfiles}}.


    A \code{\link[Biostrings]{DNAStringSet}} of the organisme of interest.
    The base pair frequency will automatically be computed and updated
    in \code{\link{genomicProfiles}} (Prefered method).
    }
}


\details{
Default value is c(0.25,0.25,0.25,0.25)
When generating a Postion Weight Matrix from a Position Frequency Matrix,
the probability of occurrence of each base pair (Base Pair Frequency) is
necessary (as originally described by Gary Stormo). It is possible to
set custom values for \code{BPFrequency} with a vector of length 4
containing the probability of occurrence of each base pair (A,C,G,T) in order.
If Base pair frequency is unknown, \code{BPFrequency} will compute base pair
frequency from a DNA sequence when building a
\code{\link{genomicProfiles}} object.
The nature of this sequence can be a\code{BSgenome} object or a
\code{\link[Biostrings]{DNAStringSet}}. In order to decrease run time,
it is advised to use \code{\link[Biostrings]{DNAStringSet}}.
}
\value{
Returns a \code{\link{genomicProfiles}} object with an updated
value for \code{BPFrequency}.
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
GPP <- genomicProfiles(PFM=PFM,PFMFormat = "JASPAR", BPFrequency=DNASequenceSet)
# Updating BPFrequency
## !! Note!! BPFrequency is used to compute PWM from PFM
## IF updated after building GPP, then it will not influence PWM
## Advised to build with BPFrequency directly
BPFrequency(GPP) <- DNASequenceSet
BPFrequency(GPP) <- c(0.25,0.25,0.25,0.25)
}
}
