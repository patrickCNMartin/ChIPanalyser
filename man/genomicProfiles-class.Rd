\name{genomicProfiles-class}
\Rdversion{1.1}
\docType{class}
\alias{genomicProfiles-class}



\title{Class \code{"genomicProfiles"}}
\description{
\code{genomicProfiles} is an S4 object serving two purposes: (i) storing internal computed data and (ii) storing paramter options. This object is parsed through the different steps of the pipeline
to facilitate that parsing and changing of paramters.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{genomicProfiles(ploidy, boundMolecules, backgroundSignal, maxSignal, lociWidth, chipMean, chipSd, chipSmooth, stepSize, noiseFilter, removeBackground, lambdaPWM, PWMpseudocount, naturalLog, noOfSites, PWMThreshold, strandRule, whichstrand, PFM, PWM, PFMFormat, BPFrequency, minPWMScore, maxPWMScore, profiles, DNASequenceLength, averageExpPWMScore)}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
  \describe{
  \item{\code{PWM}:}{Object of class \code{"matrix"}:
  A Position Weight Matrix (either supplied or internally computed if
  PFM is provided)}
  \item{\code{PFM}:}{Object of class \code{"matrix"}:
  A Position Frequency Matrix (may also be a path to file containing PFM) }
  \item{\code{PFMFormat}:}{Object of class \code{"character"}:
  A character string of one of the following: raw, transfac,JASPAR or
  sequences}
  \item{\code{BPFrequency}:}{Object of class \code{"vector"}:
  Base Pair Frequency in the genome (if a DNA sequence is provided
  (as a \code{\link[Biostrings]{DNAStringSet}} or \code{BSgenome}),
  will be automatically computed internally). Default:c(0.25,0.25,0.25,0.25)}
  \item{\code{minPWMScore}:}{Object of class \code{"vector"}:
  Lowest PWM score accros the genome (computed and updated internally) }
  \item{\code{maxPWMScore}:}{Object of class \code{"vector"}:
  Highest PWM score across the genome (computed and updated internally) }
    \item{\code{profiles}:}{Object of class \code{"GRList"}:
    Containins GRanges with sites above threshold and associated metrics
    (PWMscore and Occupancy) - Computed Internally}
    \item{\code{DNASequenceLength}:}{Object of class \code{"vector"}:
    Length of the Genome (or accesible genome) - computed internally }
    \item{\code{averageExpPWMScore}:}{Object of class \code{"vector"}:
    Average exponential PWM score across the genome
    (or accesible genome) - computed internally }
    \item{\code{ZeroBackground}:}{Object of class \code{"vector"}:
    Internal background value (computed internally)}
    \item{\code{drop}:}{Object of class \code{"vector"}:
    Stores Loci that do contain accesible DNA if it were to be the case
    (computed and updated internally) }
    \item{\code{tags}:}{Object of class \code{"character"} ~Internal Tags~ }
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
    \item{\code{lociWidth}:}{Object of class \code{"numeric"} ~~ }
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
    \item{\code{noiseFilter}:}{Object of class \code{"character"} ~Describes the
    noiseFilter method that will be applied to ChIP data (Zero, mean, median,
    sigmoid)~ }
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
\section{Extends}{
Class \code{"\linkS4class{genomicProfilesInternal}"}, directly.
Class \code{"\linkS4class{parameterOptions}"}, directly.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "genomicProfiles")}: ... }
    \item{show}{\code{signature(object = "genomicProfiles")}: ... }
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
\code{\link{genomicProfiles}}
\code{\link{parameterOptions}}
}
\examples{
showClass("genomicProfiles")
}
\keyword{classes}
