\name{parameterOptions}
\alias{parameterOptions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
parameter Options object
}
\description{
\code{parameterOptions} is an object used to store and parse the various parameters
needed throughout this analysis pipeline.
}
\usage{
parameterOptions(ploidy = 2, boundMolecules = 1000, backgroundSignal = 0, maxSignal = 1, lociWidth = 20000, chipMean = 200, chipSd = 200, chipSmooth = 250, stepSize = 10, removeBackground = 0, noiseFilter = "zero", naturalLog = TRUE, noOfSites = "all", PWMThreshold = 0.7, strandRule = "max", whichstrand = "+-", PWMpseudocount = 1, lambdaPWM = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{ploidy}{ \code{\link{ploidy}}:
A numeric Value descibing the ploidy of the organism. Default: 2}
\item{boundMolecules}{ \code{\link{boundMolecules}}:
A vector (or single value) containing the number of bound Molecules
(bound Transcription Factors): Default: 1000}
\item{backgroundSignal}{ \code{\link{backgroundSignal}}:
A numeric value descibing the ChIP-seq background Signal
(average signal from real ChIP seq data). Default: 0}
\item{maxSignal}{ \code{\link{maxSignal}}:
A numeric value describing the highest ChIP-seq signal
(from real ChIP-seq data). Default: 1}
\item{lociWidth}{ \code{\link{lociWidth}}:
A numeric value describing the width of the bins used to split ChIP profiles
parsed to \code{\link{processingChIP}}. Default = 20000}
\item{chipMean}{ \code{\link{chipMean}}:
A numeric value describing the mean width of a ChIP- seq peak: Default:200}
\item{chipSd}{ \code{\link{chipSd}}:
A numeric value describing the standard deviation of ChIP-seq peaks.
Default: 200 }
\item{chipSmooth}{ \code{\link{chipSmooth}}:
A numeric value describing the width of the window used to smooth
Occupancy profiles into ChIP profiles. Default:250}
\item{stepSize}{ \code{\link{stepSize}}:
A numeric value describing the step Size (in base pairs)
between each ChIP-seq score. Default:10 (Scored every 10 base pairs)}
\item{removeBackground}{ \code{\link{removeBackground}}:
A numeric value describing the value at which score should be removed.
Defualt:0 (If negative scores then remove) }
  \item{noiseFilter}{\code{\link{noiseFilter}}:
  A character string of one of the following: Zero, Mean, Median, or Sigmoid.
  Noise filter that will be applied to the ChIP Score during the \code{processingChIP} step.

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

\item{PWMpseudocount}{ \code{\link{PWMpseudocount}}:
A numeric value describing a PWMpseudocount for PWM computation. Default:1
}
\item{lambdaPWM}{
A vector (or single value) contaning values for the ScalingFactorPWM
(Also known as lambda).Default:1
}
}
\details{
ChIPanalyser requires a lot of parameters. \code{parameterOptions} was created with the intent
of storing and parsing these numerous arguments to the different functions. All parameters in this object are optional although strongly recommend. Some parameters are extracted and updated from
function along the pipeline e.g. maxSignal and backgroundSignal are extracted during the
\code{processingChIP} step. These paramters will be automatically parsed. If you do not which to use them ( or any other parameter) simply parse a new parameterOptions object with your desired
paramters.
}
\value{
Returns a \code{parameterOptions} with updated values.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.

}
\author{
Patrick C. N. Martin  <p.martin@essex.ac.uk>
}

\seealso{
\code{\linkS4class{genomicProfiles}}

}
\examples{
#
parameterOptions(ploidy = 2, boundMolecules = 1000, backgroundSignal = 0,
    maxSignal = 1, lociWidth = 20000, chipMean = 200, chipSd = 200,
    chipSmooth = 250, stepSize = 10, removeBackground = 0, noiseFilter = "zero",
    naturalLog = TRUE, noOfSites = "all", PWMThreshold = 0.7,
    strandRule = "max", whichstrand = "+-", PWMpseudocount = 1,
    lambdaPWM = 1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
