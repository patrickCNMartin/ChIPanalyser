\name{genomicProfilesInternal-class}
\Rdversion{1.1}
\docType{class}
\alias{genomicProfilesInternal-class}
\alias{.averageExpPWMScore<-,genomicProfilesInternal,numeric-method}
\alias{.DNASequenceLength<-,genomicProfilesInternal,vector-method}
\alias{.drop<-,genomicProfilesInternal,vector-method}
\alias{.generatePWM,genomicProfilesInternal-method}
\alias{.maxPWMScore<-,genomicProfilesInternal,vector-method}
\alias{.minPWMScore<-,genomicProfilesInternal,vector-method}
\alias{.profiles<-,genomicProfilesInternal,GRList-method}
\alias{.tags,genomicProfilesInternal-method}
\alias{.tags<-,genomicProfilesInternal,character-method}
\alias{averageExpPWMScore,genomicProfilesInternal-method}
\alias{BPFrequency,genomicProfilesInternal-method}
\alias{BPFrequency<-,genomicProfilesInternal,DNAStringSet-method}
\alias{BPFrequency<-,genomicProfilesInternal,vector-method}
\alias{DNASequenceLength,genomicProfilesInternal-method}
\alias{drop,genomicProfilesInternal-method}
\alias{maxPWMScore,genomicProfilesInternal-method}
\alias{minPWMScore,genomicProfilesInternal-method}
\alias{PFMFormat,genomicProfilesInternal-method}
\alias{PFMFormat<-,genomicProfilesInternal,character-method}
\alias{PositionFrequencyMatrix,genomicProfilesInternal-method}
\alias{PositionFrequencyMatrix<-,genomicProfilesInternal,character-method}
\alias{PositionFrequencyMatrix<-,genomicProfilesInternal,matrix-method}
\alias{PositionWeightMatrix,genomicProfilesInternal-method}
\alias{PositionWeightMatrix<-,genomicProfilesInternal,matrix-method}
\alias{profiles,genomicProfilesInternal-method}

\title{Class \code{"genomicProfilesInternal"}}
\description{
Non exported class. Represents the stripped down version of genomicProfiles.
}
\section{Objects from the Class}{
Created Internally.
}
\section{Slots}{
  \describe{
    \item{\code{PWM}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{PFM}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{PFMFormat}:}{Object of class \code{"character"} ~~ }
    \item{\code{BPFrequency}:}{Object of class \code{"vector"} ~~ }
    \item{\code{minPWMScore}:}{Object of class \code{"vector"} ~~ }
    \item{\code{maxPWMScore}:}{Object of class \code{"vector"} ~~ }
    \item{\code{profiles}:}{Object of class \code{"GRList"} ~~ }
    \item{\code{DNASequenceLength}:}{Object of class \code{"vector"} ~~ }
    \item{\code{averageExpPWMScore}:}{Object of class \code{"vector"} ~~ }
    \item{\code{ZeroBackground}:}{Object of class \code{"vector"} ~~ }
    \item{\code{drop}:}{Object of class \code{"vector"} ~~ }
    \item{\code{tags}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{.averageExpPWMScore<-}{\code{signature(object = "genomicProfilesInternal", value = "numeric")}: ... }
    \item{.DNASequenceLength<-}{\code{signature(object = "genomicProfilesInternal", value = "vector")}: ... }
    \item{.drop<-}{\code{signature(object = "genomicProfilesInternal", value = "vector")}: ... }
    \item{.generatePWM}{\code{signature(object = "genomicProfilesInternal")}: ... }
    \item{.maxPWMScore<-}{\code{signature(object = "genomicProfilesInternal", value = "vector")}: ... }
    \item{.minPWMScore<-}{\code{signature(object = "genomicProfilesInternal", value = "vector")}: ... }
    \item{.profiles<-}{\code{signature(object = "genomicProfilesInternal", value = "GRList")}: ... }
    \item{.tags}{\code{signature(object = "genomicProfilesInternal")}: ... }
    \item{.tags<-}{\code{signature(object = "genomicProfilesInternal", value = "character")}: ... }
    \item{averageExpPWMScore}{\code{signature(object = "genomicProfilesInternal")}: ... }
    \item{BPFrequency}{\code{signature(object = "genomicProfilesInternal")}: ... }
    \item{BPFrequency<-}{\code{signature(object = "genomicProfilesInternal", value = "DNAStringSet")}: ... }
    \item{BPFrequency<-}{\code{signature(object = "genomicProfilesInternal", value = "vector")}: ... }
    \item{DNASequenceLength}{\code{signature(object = "genomicProfilesInternal")}: ... }
    \item{drop}{\code{signature(object = "genomicProfilesInternal")}: ... }
    \item{maxPWMScore}{\code{signature(object = "genomicProfilesInternal")}: ... }
    \item{minPWMScore}{\code{signature(object = "genomicProfilesInternal")}: ... }
    \item{PFMFormat}{\code{signature(object = "genomicProfilesInternal")}: ... }
    \item{PFMFormat<-}{\code{signature(object = "genomicProfilesInternal", value = "character")}: ... }
    \item{PositionFrequencyMatrix}{\code{signature(object = "genomicProfilesInternal")}: ... }
    \item{PositionFrequencyMatrix<-}{\code{signature(object = "genomicProfilesInternal", value = "character")}: ... }
    \item{PositionFrequencyMatrix<-}{\code{signature(object = "genomicProfilesInternal", value = "matrix")}: ... }
    \item{PositionWeightMatrix}{\code{signature(object = "genomicProfilesInternal")}: ... }
    \item{PositionWeightMatrix<-}{\code{signature(object = "genomicProfilesInternal", value = "matrix")}: ... }
    \item{profiles}{\code{signature(object = "genomicProfilesInternal")}: ... }
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
\code{\link{genomicProfiles}}

\code{\link{parameterOptions}}

}
\examples{
showClass("genomicProfilesInternal")
}
\keyword{classes}
