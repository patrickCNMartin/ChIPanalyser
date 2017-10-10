\name{stepSize<-}
\alias{stepSize<-}
\alias{stepSize<--methods}

\title{
Setter Method for the \code{stepSize} slot in a
\code{\link{occupancyProfileParameters}}
}
\description{
Setter Method for the \code{stepSize} slot in a
\code{\link{occupancyProfileParameters}}
}
\usage{
stepSize(object) <- value
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{occupancyProfileParameters}} object
    }
    \item{value}{
    \code{value} is a positive numeric value that will be assigned to the
    \code{stepSize} slot in a \code{\link{occupancyProfileParameters}} object.
    Default is set at 10 base pairs.
    }
}
\details{
It possible to restrict the size of the ChIP-seq-like profile produced by
\code{\link{computeChipProfile}}. Instead of returning ChIP-seq like score
for each base pair, it is possible to skip base pairs and only return the
predicted enrichement score for every "n" base pair
(n is the value assigned to stepSize). This will reduce the size of the
output data (unless step size is very large, this will not affect the
accuracy of the model). Default is set at 10 base pairs.
}
\value{
Returns a \code{occupancyProfileParameters} object with an updated value
for the \code{stepSize} slot.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>

}

\examples{
# Building occupancyProfileParameters object
OPP <- occupancyProfileParameters()
#Setting new Value for maxSignal
stepSize(OPP) <- 20

}
