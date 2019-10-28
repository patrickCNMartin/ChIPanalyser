\name{processingChIP}
\alias{processingChIP}

\title{
Pre-processing ChIP-seq data
}
\description{
\code{processingChIP}
will process and extract ChIP scores at a set of loci of interest.
}

\usage{
processingChIP(profile,loci=NULL,reduce=NULL,
     peaks=NULL,chromatinState=NULL,parameterOptions=NULL,
     cores=1)
}

\arguments{
    \item{profile}{
    \code{profile} is a path to a UCSC format file, a GRanges or data frame. The input data  should contain 4 columns:chromosome, start , end and score.
    }
    \item{loci}{
    \code{loci} is \code{\link[GenomicRanges]{GRanges}} describing the loci
    at which ChIP scores should be extracted.
    If NULL, a set of Loci will extracted from \code{profile}. The data provided will then be split into bins of width equal to \code{lociWidth} (Default 20kbp)
     Default=NULL
    }
    \item{reduce}{
    \code{reduce} is a the top regions to select based on the mean ChIP score.
    If \code{peaks} are provided, regions overlappling with known peaks will be
    selected based on highest ChIP score. If NULL, all regions will
    be considered. Default=NULL

    }
    \item{parameterOptions}{
    \code{parameterOptions} is an
    \code{\link{parameterOptions}} object containing chip
    Parameters to be parsed for ChIP score extraction. If NULL,
    \code{\link{parameterOptions}} will be built internally with
    default ChIP extraction parameters (see  \code{\link{chipSmooth}},
    \code{\link{chipSd}} and \code{\link{chipMean}})  Default=NULL
    }
    \item{peaks}{
    \code{peaks} is a path to UCSC format file or a GRanges object containing
    location of ChIP peaks.  Default=NULL
    }
    \item{chromatinState}{
    \code{chromatinState} is a \code{\link[GenomicRanges]{GRanges}} containing
    Accessible DNA or chromatin States If provided, regions will be selected only if they
    contain accesible DNA.  Default=NULL
    }

    \item{cores}{
    \code{cores} is the number of cores used to extract ChIP scores. Default = 1
    }
}
\details{
When using \code{\link{computeOptimal}}, it is required to supply real ChIP
data in order to have a point of comparison. The corralation and MSE Scores are
computed based of how well the model fits biological data.
\code{processingChIP} will extract this data from ChIP data at loci
of interest. When using the \code{reduce} option, this function will only
select the top regions based on peak height or mean ChIP score.
\code{processingChIP} will also extract maxSignal and backgroundSignal from
ChIP data and parse it to an \code{\link{parameterOptions}} object.
}
\value{
Returns a ChIPScore object containing extracted (and normalised) ChIP scores, the loci
of interest and newly extracted Parameters(e.g. maxSignal)
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

ChIP<-processingChIP(profile=eveLocusChip, loci=eveLocus)

}
