\name{computeOptimal}
\alias{computeOptimal}

\title{
compute Optimal Parameters
}
\description{
\code{ChIPanalyser} contains a set of functions some of which require two
parameters known as \code{\link{lambdaPWM}} and as
\code{\link{boundMolecules}}. These two paramters are not always known.
\code{computeOptimal} will compute these values by maximising the
correlation and minimising the Mean Squared Error between a predicted
ChIP-seq-like profile and a real ChIP-seq profile for a given loci.
}
\usage{
computeOptimal(genomicProfiles,DNASequenceSet, ChIPScore,chromatinState = NULL,
    parameterOptions = NULL, optimalMethod = "all",rank=FALSE,returnAll=TRUE,
    peakMethod="moving_kernel",cores=1)
}

\arguments{
     \item{genomicProfiles}{
     \code{genomicProfiles} is a \code{\link{genomicProfiles}}
     object containing at least a Postion Frequency Matrix or a
     Position Weight Matrix. It is strongly advised to customize this object to
     increase goodness of fit of the model when compared to real ChIP-seq data.
     }
    \item{DNASequenceSet}{
    \code{DNASequenceSet} is a \code{\link[Biostrings]{DNAStringSet}} or
    a \code{BSgenome} of the full sequence of the organism of interest.
    }

    \item{ChIPScore}{
    \code{ChIPScore} is a named list containing ChIP-seq enrichements for
    each Loci of interest. This Profile should be normalised to a base
    pair level. In other words, there should be an enrichement score for
    each base pair of a given Locus.
    }


    \item{chromatinState}{
    \code{chromatinState} is a \code{\link[GenomicRanges]{GRanges}} object
    containing either accesible sites or DNA affinity scores.
    }

    \item{parameterOptions}{
    \code{\link{parameterOptions}} is a
    \code{\link{parameterOptions}} object.
    If this object is not provided (parameterOptions = NULL),
    a new object will be created internally.
    However, it is strongly advised to tailor this object to maximise the
    goodness of fit of the model when compared to ChIP-seq data.
    }
    \item{optimalMethod}{
    \code{optimalMethod} is a character string which determines which
    method for optimal parameter selection should be selected. optimalMethod
    can be one of the following:
    pearson, spearman, kendall, ks, fscore, geometric,MSE, or all. Default is set
    to all.
    }
    \item{rank}{
    \code{rank} is a logical value indicating if optimal parameters should be
    based on rank (parameter combination occuring the most over all regions) or
    avaerage score (best perfomring combination of paramters on average over all
    regions selected). DEFAULT = FALSE
    }
    \item{returnAll}{
    \code{returnAll} is a logical value indicating if all internal objects
    should be returned. DEFAULT = TRUE. Internal objects are the following:
    Occupancy Scores, ChIP like profiles, goodness of fit metrics and optimal
    paramters. If set to FALSE, \code{computeOptimal} will only return the optimal
    parameters.
    }
    \item{peakMethod}{
    \code{peakMethod} is a character string  of one of the following:
    c("moving_kernel","truncated_kernel","exact"). If set to moving_kernel,
    the peaks will be approximated using Rcpp (Default).
    If set to truncated_kernel, the peaks will be approximated
    however this method does not require Rcpp. If set to exact, the peaks will
    not be approximated.

    }
    \item{cores}{
    \code{cores} is the number cores that will be used to compute optimal
    set of parameters.
    }
}
\details{
In order to backward infer the values of \code{\link{lambdaPWM}}
and \code{\link{boundMolecules}}, it is possible to use the
\code{computeOptimal} to find these parameters.
It should be noted that this functions requires a ChIP-seq data input.
\code{ChIPScore} (ChIP-seq data). This should be the output of the \code{\link{processingChIP}} function.
}
\value{
\code{computeOptimal} returns a list respectivly described as the optimal
set of Parameters (lambda - \code{\link{lambdaPWM}} and
\code{\link{boundMolecules}}), the optimal matrix (a matrix containing
accuracy estimates dependant on the parameter chosen), and finally the
chosen parameter. If the parameter that was chosen was "all",
then each element of this list will contain the optimal set of
parameters, optimal matricies for all of the aforementioned paramters (see \code{optimalMethod}).
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.

}
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
eveLocusChip<-processingChIP(eveLocusChip,eveLocus)
#Building data objects
GPP <- genomicProfiles(PFM=PFM,BPFrequency=DNASequenceSet)
OPP <- parameterOptions()
#Computing Optimal set of Parameters
optimalParam <- computeOptimal(genomicProfiles = GPP,
    DNASequenceSet = DNASequenceSet,
    ChIPScore = eveLocusChip,
    chromatinState = Access,
    parameterOptions = OPP,
    parameter = "all",
    peakMethod="moving_kernel")
}
}
