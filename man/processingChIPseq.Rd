\name{processingChIPseq}
\alias{processingChIPseq}

\title{
Pre-processing ChIP-seq data from UCSC format file
}
\description{
\code{processingChIPseq}
will process and extract ChIP scores at a set of loci of interest.
}

\usage{
processingChIPseq(profile,loci=NULL,reduce=NULL,
occupancyProfileParameters=NULL,
peaks=NULL,Access=NULL,noiseFilter=c("zero","mean","median","sigmoid"),cores=1)
}

\arguments{
    \item{profile}{
    \code{profile} is a path to a UCSC format file, a GRanges or a data.frame
    containing ChIP scores. The input data frame should contain 4 columns:
    chromosome, start , end and score. This is also applicable for the GRanges
    format.
    }
    \item{loci}{
    \code{loci} is \code{\link[GenomicRanges]{GRanges}} describing the loci
    at which ChIP scores should be extracted.
    If NULL, a set of Loci will extracted from \code{profile} based on
    chromosomes. However, we STRONGLY recommend to use a
    \code{\link[GenomicRanges]{GRanges}} of loci of interest. Default=NULL
    }
    \item{reduce}{
    \code{reduce} is a the top regions to select based on the mean ChIP score.
    If \code{peaks} are provided, regions overlappling with known peaks will be
    selected based on highest ChIP score. If NULL, all regions will
    be considered. Default=NULL

    }
    \item{occupancyProfileParameters}{
    \code{occupancyProfileParameters} is an
    \code{\link{occupancyProfileParameters}} object containing chip
    Parameters to be parsed for ChIP score extraction. If NULL,
    \code{\link{occupancyProfileParameters}} will be built internally with
    default ChIP extraction parameters (see  \code{\link{chipSmooth}},
    \code{\link{chipSd}} and \code{\link{chipMean}})  Default=NULL
    }
    \item{peaks}{
    \code{peaks} is a path to UCSC format file or a GRanges object containing
    location of ChIP peaks.  Default=NULL
    }
    \item{Access}{
    \code{Access} is a \code{\link[GenomicRanges]{GRanges}} containing
    Accessible DNA. If provided, regions will be selected only if they
    contain accesible DNA.  Default=NULL
    }
    \item{noiseFilter}{\code{noiseFilter}: Noise filtering method that
    should be used on ChIP-seq data. Four methods are available: Zero, Mean,
    Median and Sigmoid. Zero removes all ChIP-seq scores bellow zero, mean under
    the mean score, median under median score and sigmoid assignes a weight
    to each score based on a logistic regression curve. Mid point is set at 95
    95 quantile of ChIP-seq scores. Below midpoint will receive a score between
    0 and 1 , everything above will receive a score between 1 and 2}
    \item{cores}{
    \code{cores} is the number of cores used to extract ChIP scores. Default = 1
    }
}
\details{
When using \code{\link{computeOptimal}}, it is required to supply real ChIP
data in order to have a point of comparison. The corralation and MSE Scores are
computed based of how well the model fits biological data.
\code{processingChIPseq} will extract this data from ChIP data at loci
of interest. When using the \code{reduce} option, this function will only
select the top regions based on peak height or mean ChIP score.
\code{processingChIPseq} will also extract maxSignal and backgroundSignal from
ChIP data and parse it to an \code{\link{occupancyProfileParameters}} object.
}
\value{
If using reduce, will return a list of two elements.
The first element will contain a list with extracted ChIP data and a new set
of top scoring loci loci. The second element will contain a
\code{\link{occupancyProfileParameters}} object with maxSignal
and backgroundSignal slot updated.
If NOT using reduce, the first element will only contain
ChIP score at loci of interest the second will still contain  a
\code{\link{occupancyProfileParameters}} object with maxSignal
and backgroundSignal slot updated.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}
\author{
Patrick C.N. Martin <pm16057@essex.ac.uk>
}
\examples{


#Data extraction
data(ChIPanalyserData)

## Extracting ChIP scores at loci of interest

ChIP<-processingChIPseq(profile=eveLocusChip, loci=eveLocus)

}
