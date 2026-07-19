\name{profiles-methods}
\docType{methods}
\alias{profiles}
\title{ ~~ Methods for Function \code{profiles}  ~~}
\description{
Accessor method for \code{profiles} in a \code{\link{genomicProfiles}} object
}
\value{
 A \code{GRanges}, \code{GRangesList}, or \code{list} containing the computed PWM score, Occupancy, or ChIP-like profiles.
}

\section{Methods}{
\describe{

\item{\code{profiles(oject)}}{
Computed PWM scores, Occupancy or ChIP-seq like profiles for loci of interest and
paramter combination of interest.
}
}}
\examples{
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BEAF-32.pfm")
GPP <- genomicProfiles(PFM=PFM,PFMFormat="JASPAR")
profiles(GPP)
}
\keyword{methods}
