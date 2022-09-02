\name{parameterOptions-class}
\Rdversion{1.1}
\docType{class}
\alias{parameterOptions-class}
\alias{.paramTag,parameterOptions-method}
\alias{.paramTag<-,parameterOptions,character-method}
\alias{.ZeroBackground,parameterOptions-method}
\alias{.ZeroBackground<-,parameterOptions,vector-method}
\alias{backgroundSignal,parameterOptions-method}
\alias{backgroundSignal<-,parameterOptions,numeric-method}
\alias{boundMolecules,parameterOptions-method}
\alias{boundMolecules<-,parameterOptions,vector-method}
\alias{chipMean,parameterOptions-method}
\alias{chipMean<-,parameterOptions,numeric-method}
\alias{chipSd,parameterOptions-method}
\alias{chipSd<-,parameterOptions,numeric-method}
\alias{chipSmooth,parameterOptions-method}
\alias{chipSmooth<-,parameterOptions,vector-method}
\alias{lambdaPWM,parameterOptions-method}
\alias{lambdaPWM<-,parameterOptions,vector-method}
\alias{lociWidth,parameterOptions-method}
\alias{lociWidth<-,parameterOptions,numeric-method}
\alias{maxSignal,parameterOptions-method}
\alias{maxSignal<-,parameterOptions,numeric-method}
\alias{naturalLog,parameterOptions-method}
\alias{naturalLog<-,parameterOptions,logical-method}
\alias{noiseFilter,parameterOptions-method}
\alias{noiseFilter<-,parameterOptions,character-method}
\alias{noOfSites,parameterOptions-method}
\alias{noOfSites<-,parameterOptions,character-method}
\alias{noOfSites<-,parameterOptions,numeric-method}
\alias{ploidy,parameterOptions-method}
\alias{ploidy<-,parameterOptions,numeric-method}
\alias{PWMpseudocount,parameterOptions-method}
\alias{PWMpseudocount<-,parameterOptions,numeric-method}
\alias{PWMThreshold,parameterOptions-method}
\alias{PWMThreshold<-,parameterOptions,numeric-method}
\alias{removeBackground,parameterOptions-method}
\alias{removeBackground<-,parameterOptions,vector-method}

\alias{stepSize,parameterOptions-method}
\alias{stepSize<-,parameterOptions,numeric-method}
\alias{strandRule,parameterOptions-method}
\alias{strandRule<-,parameterOptions,character-method}
\alias{whichstrand,parameterOptions-method}
\alias{whichstrand<-,parameterOptions,character-method}

\title{Class \code{"parameterOptions"}}
\description{
\code{parameterOptions} is an object used to store and parse the various parameters
needed throughout this analysis pipeline.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{parameterOptions(ploidy, boundMolecules, backgroundSignal, maxSignal, lociWidth, chipMean, chipSd, chipSmooth, stepSize, noiseFilter, removeBackground, lambdaPWM, PWMpseudocount, naturalLog, noOfSites, PWMThreshold, strandRule, whichstrand)}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
  \describe{
  \item{\code{ploidy}:}{Object of class \code{"numeric"}:
  A numeric Value descibing the ploidy of the organism. Default: 2 }
  \item{\code{boundMolecules}:}{Object of class \code{"vector"}:
  A vector (or single value) containing the number of bound Molecules
  (bound Transcription Factors): Default: 1000 }
  \item{\code{backgroundSignal}:}{Object of class \code{"numeric"}:
  A numeric value descibing the ChIP-seq background Signal
  (average signal from real ChIP seq data). Default: 0}
  \item{\code{maxSignal}:}{Object of class \code{"numeric"}:
  A numeric value describing the highest ChIP-seq signal
  (from real ChIP-seq data). Default: 1}
  \item{\code{lociWidth}:}{Object of class \code{"numeric"}:
  A numeric value describing bin size when splitting ChIP seq scores). Default: 20 000}
  \item{\code{chipMean}:}{Object of class \code{"numeric"}:
  A numeric value describing the mean width of a ChIP- seq peak.
  Default:150 }
  \item{\code{chipSd}:}{Object of class \code{"numeric"}:
  A numeric value describing the standard deviation of ChIP-seq peaks.
  Default: 150 }
  \item{\code{chipSmooth}:}{Object of class \code{"vector"}:
  A numeric value describing the width of the window used to smooth
  Occupancy profiles into ChIP profiles. Default:250 }
  \item{\code{stepSize}:}{Object of class \code{"numeric"}:
  A numeric value describing the step Size (in base pairs) between
  each ChIP-seq score. Default:10 (Scored every 10 base pairs)}
  \item{\code{removeBackground}:}{Object of class \code{"numeric"}:
  A numeric value describing the value at which score should be removed.
  Defualt:0 (If negative scores then remove) }
    \item{\code{noiseFilter}:}{Object of class \code{"character"} Describes
    noiseFilter method applied to ChIP scores}
    \item{\code{PWMThreshold}:}{Object of class \code{"numeric"}:
    Threshold at which PWM Score should be selected (only sites above
    threshold will be selected - between 0 and 1) }
    \item{\code{strandRule}:}{Object of class \code{"character"}:
    "mean", "max" or "sum"  will dertermine how strand should be handle
    for computing PWM Scores. Default : "max" }
    \item{\code{whichstrand}:}{Object of class \code{"character"}:
    "+","-" or "+-" on which strand should PWM Score be computed.
    Default: "+-" }
    \item{\code{lambdaPWM}:}{Object of class \code{"vector"} A vector (or single value) contaning values for lambdaPWM Default:1}
    \item{\code{naturalLog}:}{Object of class \code{"logical"}:
    A logical value describing if natural Log will be used to compute
    the PWM (if FALSE then log2 will be used). Default: TRUE }
    \item{\code{noOfSites}:}{Object of class \code{"nos"} A Positive integer descibing number of sites (in base pair) should
    be used from the PFM to compute PWM. Default =0 (Full width of
    binding site will be used when set to 0) }
    \item{\code{PWMpseudocount}:}{Object of class \code{"numeric"}:
    A numeric value describing a PWMpseudocount for PWM computation.
    Default:1 }
    \item{\code{paramTag}:}{Object of class \code{"character"} ~Internal~ }
  }
}
\section{Methods}{
  \describe{
    \item{.paramTag}{\code{signature(object = "parameterOptions")}: ... }
    \item{.paramTag<-}{\code{signature(object = "parameterOptions", value = "character")}: ... }
    \item{.ZeroBackground}{\code{signature(object = "parameterOptions")}: ... }
    \item{.ZeroBackground<-}{\code{signature(object = "parameterOptions", value = "vector")}: ... }
    \item{backgroundSignal}{\code{signature(object = "parameterOptions")}: ... }
    \item{backgroundSignal<-}{\code{signature(object = "parameterOptions", value = "numeric")}: ... }
    \item{boundMolecules}{\code{signature(object = "parameterOptions")}: ... }
    \item{boundMolecules<-}{\code{signature(object = "parameterOptions", value = "vector")}: ... }
    \item{chipMean}{\code{signature(object = "parameterOptions")}: ... }
    \item{chipMean<-}{\code{signature(object = "parameterOptions", value = "numeric")}: ... }
    \item{chipSd}{\code{signature(object = "parameterOptions")}: ... }
    \item{chipSd<-}{\code{signature(object = "parameterOptions", value = "numeric")}: ... }
    \item{chipSmooth}{\code{signature(object = "parameterOptions")}: ... }
    \item{chipSmooth<-}{\code{signature(object = "parameterOptions", value = "vector")}: ... }
    \item{initialize}{\code{signature(.Object = "parameterOptions")}: ... }
    \item{lambdaPWM}{\code{signature(object = "parameterOptions")}: ... }
    \item{lambdaPWM<-}{\code{signature(object = "parameterOptions", value = "vector")}: ... }
    \item{lociWidth}{\code{signature(object = "parameterOptions")}: ... }
    \item{lociWidth<-}{\code{signature(object = "parameterOptions", value = "numeric")}: ... }
    \item{maxSignal}{\code{signature(object = "parameterOptions")}: ... }
    \item{maxSignal<-}{\code{signature(object = "parameterOptions", value = "numeric")}: ... }
    \item{naturalLog}{\code{signature(object = "parameterOptions")}: ... }
    \item{naturalLog<-}{\code{signature(object = "parameterOptions", value = "logical")}: ... }
    \item{noiseFilter}{\code{signature(object = "parameterOptions")}: ... }
    \item{noiseFilter<-}{\code{signature(object = "parameterOptions", value = "character")}: ... }
    \item{noOfSites}{\code{signature(object = "parameterOptions")}: ... }
    \item{noOfSites<-}{\code{signature(object = "parameterOptions", value = "character")}: ... }
    \item{noOfSites<-}{\code{signature(object = "parameterOptions", value = "numeric")}: ... }
    \item{ploidy}{\code{signature(object = "parameterOptions")}: ... }
    \item{ploidy<-}{\code{signature(object = "parameterOptions", value = "numeric")}: ... }
    \item{PWMpseudocount}{\code{signature(object = "parameterOptions")}: ... }
    \item{PWMpseudocount<-}{\code{signature(object = "parameterOptions", value = "numeric")}: ... }
    \item{PWMThreshold}{\code{signature(object = "parameterOptions")}: ... }
    \item{PWMThreshold<-}{\code{signature(object = "parameterOptions", value = "numeric")}: ... }
    \item{removeBackground}{\code{signature(object = "parameterOptions")}: ... }
    \item{removeBackground<-}{\code{signature(object = "parameterOptions", value = "vector")}: ... }
    \item{show}{\code{signature(object = "parameterOptions")}: ... }
    \item{stepSize}{\code{signature(object = "parameterOptions")}: ... }
    \item{stepSize<-}{\code{signature(object = "parameterOptions", value = "numeric")}: ... }
    \item{strandRule}{\code{signature(object = "parameterOptions")}: ... }
    \item{strandRule<-}{\code{signature(object = "parameterOptions", value = "character")}: ... }
    \item{whichstrand}{\code{signature(object = "parameterOptions")}: ... }
    \item{whichstrand<-}{\code{signature(object = "parameterOptions", value = "character")}: ... }
	 }
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
Patrick C.N. Martin and Nicolae Radu Zabe (2020) Dissecting the binding mechanisms of transcription factors to DNA using a statistical thermodynamics framework. CSBJ, 18, 3590-3605.
}
\author{
Partick C. N. Martin <p.martin@essex.ac.uk>
}

\seealso{
\code{\linkS4class{genomicProfiles}}
}
\examples{
showClass("parameterOptions")
}
\keyword{classes}
