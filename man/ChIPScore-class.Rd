\name{ChIPScore-class}
\Rdversion{1.1}
\docType{class}

\alias{ChIPScore-class}
\alias{loci,ChIPScore-method}
\alias{scores,ChIPScore-method}


\title{Class \code{"ChIPScore"}}
\description{
 \code{ChIPScore} is the result of the \code{\link{processingChIP}} function. This object contains the extracted ChIP Score from ChIP data, the loci of interest and optional paramters associated to ChIPanalyser. The loci of interest will either be user provided or the top n regions as defined by the reduce argument im \code{processingChIP}. This object has the sole purpose of aiding the storage and parsing of data and parameters.
}
\section{Objects from the Class}{
Object of this class are created internaly and will be parsed to other objects as is.
}
\section{Slots}{
  \describe{
    \item{\code{scores}:}{Object of class \code{"list"} List of extracted ChIP scores }
    \item{\code{loci}:}{Object of class \code{"loci"} GRanges containing loci of interest  }
    \item{\code{ploidy}:}{Object of class \code{"numeric"} Ploidy level of the organism }
    \item{\code{boundMolecules}:}{Object of class \code{"vector"} Number of Bound molecules to the DNA }
    \item{\code{backgroundSignal}:}{Object of class \code{"numeric"} ChIP background signal (average ChIP score) }
    \item{\code{maxSignal}:}{Object of class \code{"numeric"} max ChIP signal }
    \item{\code{lociWidth}:}{Object of class \code{"numeric"} Width of loci if reduce is used and no loci are provided }
    \item{\code{chipMean}:}{Object of class \code{"numeric"} Average ChIP peak width  }
    \item{\code{chipSd}:}{Object of class \code{"numeric"} Standard Deviation of ChIP peak width  }
    \item{\code{chipSmooth}:}{Object of class \code{"vector"} Smoothing window width for ChIP score  }
    \item{\code{stepSize}:}{Object of class \code{"numeric"} Defining resolution size of ChIP like profiles (10bp = signal will be only considered every 10bp) }
    \item{\code{removeBackground}:}{Object of class \code{"numeric"} Signal Threshold to be removed. Default removes all negative scores}
    \item{\code{noiseFilter}:}{Object of class \code{"character"} Type of noise filter to be used on ChIP data.}
    \item{\code{PWMThreshold}:}{Object of class \code{"numeric"} Threshold of PWM scores that will be selected }
    \item{\code{strandRule}:}{Object of class \code{"character"} Rule to compute strand score (max, mean or sum) }
    \item{\code{whichstrand}:}{Object of class \code{"character"} Which strand should be used to compute PWM scores.}
    \item{\code{lambdaPWM}:}{Object of class \code{"vector"} Lambda value - Scaling factor to the PWM }
    \item{\code{naturalLog}:}{Object of class \code{"logical"} PFM to PWM conversion log transform ( natural log or log2) }
    \item{\code{noOfSites}:}{Object of class \code{"nos"} Number of Sites in the PWM that should be used to compute PWM scores. }
    \item{\code{PWMpseudocount}:}{Object of class \code{"numeric"} PWM pseudocount value for PFM to PWM conversion. }
    \item{\code{paramTag}:}{Object of class \code{"character"} Internal Tag - Code progression }
  }
}
\section{Extends}{
Class \code{"\linkS4class{parameterOptions}"}, directly.
}
\section{Methods}{
  \describe{
    \item{.loci<-}{\code{signature(object = "ChIPScore", value = "loci")}: ... }
    \item{.scores<-}{\code{signature(object = "ChIPScore", value = "list")}: ... }
    \item{initialize}{\code{signature(.Object = "ChIPScore")}: ... }
    \item{loci}{\code{signature(object = "ChIPScore")}: ... }
    \item{scores}{\code{signature(object = "ChIPScore")}: ... }
    \item{show}{\code{signature(object = "ChIPScore")}: ... }
	 }
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
Patrick C.N. Martin and Nicolae Radu Zabe (2020) Dissecting the binding mechanisms of transcription factors to DNA using a statistical thermodynamics framework. CSBJ, 18, 3590-3605.
}
\author{
Patrick C.N. Martin
}



\seealso{
\code{\link{processingChIP}}
}
\examples{
showClass("ChIPScore")
}
\keyword{classes}
