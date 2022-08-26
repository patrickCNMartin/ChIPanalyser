\name{BPFrequency}

\alias{BPFrequency}


\title{
Accessor method for \code{BPFrequency} slot in a
\code{\link{genomicProfiles}} object.
}
\description{
Extract or Access \code{BPFrequency} slot in a
\code{\link{genomicProfiles}} object.
}
\usage{
BPFrequency(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfiles}}
    }
}
\details{
Default value is c(0.25,0.25,0.25,0.25)
When generating a Postion Weight Matrix from a Position Frequency Matrix,
the probability
of occurrence of each base pair (Base Pair Frequency) is necessary
(as originally described by Gary Stormo). It is possible to
set custom values for \code{BPFrequency} with a vector of length 4
containing the probability of occurrence of each base pair (A,C,G,T) in order.
If Base pair frequency is unknown, \code{BPFrequency} will compute base pair
frequency from a DNA sequence. The nature of this sequence can be a
\code{BSgenome} or a \code{\link[Biostrings]{DNAStringSet}}.
In order to decrease run time, it is advised to use
\code{\link[Biostrings]{DNAStringSet}}
}
\value{
Returns \code{BPFrequency} slot in
\code{\link{genomicProfiles}} object.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.

}
\author{
Patrick C.N. Martin <pcnmartin@gmail.com>
}
\examples{
data(ChIPanalyserData)
# path to Position Frequency Matrix
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BEAF-32.pfm")
# Building genomicProfiles object
GPP <- genomicProfiles(PFM=PFM,PFMFormat="JASPAR")
#Extracting BPFrequency slot
BPFrequency(GPP)
}
