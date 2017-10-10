\name{occupancyProfileParameters-class}
\Rdversion{1.1}
\docType{class}
\alias{occupancyProfileParameters-class}

\alias{backgroundSignal<-,occupancyProfileParameters,numeric-method}
\alias{backgroundSignal,occupancyProfileParameters-method}
\alias{boundMolecules<-,occupancyProfileParameters,vector-method}
\alias{boundMolecules,occupancyProfileParameters-method}
\alias{chipMean<-,occupancyProfileParameters,numeric-method}
\alias{chipMean,occupancyProfileParameters-method}
\alias{chipSd<-,occupancyProfileParameters,numeric-method}
\alias{chipSd,occupancyProfileParameters-method}
\alias{chipSmooth<-,occupancyProfileParameters,vector-method}
\alias{chipSmooth,occupancyProfileParameters-method}
\alias{maxSignal<-,occupancyProfileParameters,numeric-method}
\alias{maxSignal,occupancyProfileParameters-method}
\alias{ploidy<-,occupancyProfileParameters,numeric-method}
\alias{ploidy,occupancyProfileParameters-method}
\alias{removeBackground<-,occupancyProfileParameters,vector-method}
\alias{removeBackground,occupancyProfileParameters-method}
\alias{show,occupancyProfileParameters-method}
\alias{stepSize<-,occupancyProfileParameters,numeric-method}
\alias{stepSize,occupancyProfileParameters-method}
\alias{thetaThreshold<-,occupancyProfileParameters,numeric-method}
\alias{thetaThreshold,occupancyProfileParameters-method}

\title{Class \code{"occupancyProfileParameters"}}
\description{
\code{occupancyProfileParameters} is a parameter storage object.
Similarly to \code{\link{genomicProfileParameters}},
this object serves the purpose of handling the set of parameters needed
when using certain functions.
}
\section{Objects from the Class}{
Objects can be created by calls of the form
\code{occupancyProfileParameters(ploidy = 2 , boundMolecules = 1000 ,
    backgroundSignal = 0 ,  maxSignal = 1 , chipMean = 150 ,   chipSd = 150 ,
    chipSmooth = 250 , stepSize = 10 ,  removeBackground = 0 ,
    thetaThreshold = 0.1 )}

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
    \item{\code{thetaThreshold}:}{Object of class \code{"numeric"}:
    A numeric value describing a threshold at which low ratios of correlation
    over Mean Squared Error (
    for accuracy estimation) should be removed. Default: 0.1 }
    }
}
\section{Methods}{
    \describe{
    \item{backgroundSignal<-}{\code{signature(object =
    "occupancyProfileParameters", value = "numeric")}: Slot Setter Method  }
    \item{backgroundSignal}{\code{signature(object =
    "occupancyProfileParameters")}: Slot Accessor Method }
    \item{boundMolecules<-}{\code{signature(object =
    "occupancyProfileParameters", value = "vector")}: Slot Setter Method }
    \item{boundMolecules}{\code{signature(object =
    "occupancyProfileParameters")}: Slot Accessor Method }
    \item{chipMean<-}{\code{signature(object =
    "occupancyProfileParameters", value = "numeric")}: Slot Setter Method }
    \item{chipMean}{\code{signature(object =
    "occupancyProfileParameters")}: Slot Accessor Method }
    \item{chipSd<-}{\code{signature(object =
    "occupancyProfileParameters", value = "numeric")}: Slot Setter Method }
    \item{chipSd}{\code{signature(object =
    "occupancyProfileParameters")}: Slot Accessor Method }
    \item{chipSmooth<-}{\code{signature(object =
    "occupancyProfileParameters", value = "vector")}: Slot Setter Method }
    \item{chipSmooth}{\code{signature(object =
    "occupancyProfileParameters")}:Slot Accessor Method }
    \item{maxSignal<-}{\code{signature(object =
    "occupancyProfileParameters", value = "numeric")}: Slot Setter Method. }
    \item{maxSignal}{\code{signature(object =
    "occupancyProfileParameters")}: Slot Accessor Method }
    \item{ploidy<-}{\code{signature(object =
    "occupancyProfileParameters", value = "numeric")}: Slot Setter Method }
    \item{ploidy}{\code{signature(object =
    "occupancyProfileParameters")}: Slot Accessor Method }
    \item{removeBackground<-}{\code{signature(object =
    "occupancyProfileParameters", value = "vector")}:Slot Setter Method }
    \item{removeBackground}{\code{signature(object =
    "occupancyProfileParameters")}: Slot Setter Method }
    \item{show}{\code{signature(object =
    "occupancyProfileParameters")}: Slot Accessor Method }
    \item{stepSize<-}{\code{signature(object =
    "occupancyProfileParameters", value = "numeric")}: Slot Setter Method }
    \item{stepSize}{\code{signature(object =
    "occupancyProfileParameters")}: Slot Accessor Method }
    \item{thetaThreshold<-}{\code{signature(object =
    "occupancyProfileParameters", value = "numeric")}: Slot Setter Method}
    \item{thetaThreshold}{\code{signature(object =
    "occupancyProfileParameters")}: Slot Accessor Method }
    }
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>

}

\seealso{
\code{\link{occupancyProfileParameters}}
\code{\link{occupancyProfileParameters-class}}
\code{\link{genomicProfileParameters}}

}
\examples{
OPP <- occupancyProfileParameters()
# OR
OPP <- occupancyProfileParameters(ploidy=2,
    boundMolecules=5000,backgroundSignal=0)
}

\keyword{classes}
