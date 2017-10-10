\name{maxSignal}
\alias{maxSignal}
\alias{maxSignal-methods}

\title{
Accessor method for the \code{maxSignal} slot in a
\code{\link{occupancyProfileParameters}} object.
}
\description{
Accessor method for the \code{maxSignal} slot in a
\code{\link{occupancyProfileParameters}} object.
}
\usage{
maxSignal(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{occupancyProfileParameters}} object.
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
Returns the value assigned to the \code{maxSignal} slot in a
\code{\link{occupancyProfileParameters}} object.
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
maxSignal(OPP)

}
