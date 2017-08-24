\name{maxSignal<-}
\alias{maxSignal<-}
\alias{maxSignal<--methods}

\title{
Setter method for \code{maxSignal} slot in a
\code{\link{occupancyProfileParameters}} object.
}
\description{
Setter method for \code{maxSignal} slot in a
\code{\link{occupancyProfileParameters}} object.

}
\usage{
maxSignal(object) <- value
}
\arguments{
    \item{object}{
    \code{object} is a \code{\link{occupancyProfileParameters}} object.
    }
    \item{value}{
    \code{value} is a numerical value to be assigned to the
    \code{maxSignal} slot.
    }
}
\details{
In the context of \code{ChIPanalyser}, \code{maxSignal} represents the
maximum normalised ChIP-Seq signal of a given Transcription factor
(or DNA binding protein). Although, A default value of 1 has been assigned to
this slot, we strongly recommend to tailor this value accordingly.
We strongly encourage to set values when building a
\code{\link{occupancyProfileParameters}} object.
}
\value{
Returns a \code{\link{occupancyProfileParameters}} with an updated
value for \code{maxSignal}.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.

}
\author{
Patrick C.N. Martin <p.martin@essex.ac.uk>
}
\examples{
# Building occupancyProfileParameters object
OPP <- occupancyProfileParameters()
#Setting new Value for maxSignal
maxSignal(OPP) <- 1.8

}
