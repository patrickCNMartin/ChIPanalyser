\name{stepSize}


\alias{stepSize}
\title{
Accessor method of the \code{stepSize} slot in
\code{\linkS4class{parameterOptions}} object
}
\description{
Accessor method of the \code{stepSize} slot in
\code{\linkS4class{parameterOptions}} object
}
\usage{
stepSize(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\linkS4class{parameterOptions}} object.
    }
}
\details{
It possible to restrict the size of the ChIP-seq-like profile produced
by \code{\link{computeChIPProfile}}. Instead of returning ChIP-seq like
score for each base pair, it is possible to skip base pairs and only
return the predicted enrichement score for every "n" base pair
(n is the value assigned to stepSize). This will reduce the size of the
output data (unless step size is very large, this will not affect
the accuracy of the model). Default is set at 10 base pairs.

}
\value{
Returns the value assigned to the \code{stepSize} slot in a
\code{\linkS4class{parameterOptions}}
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>

}
\examples{
# Building parameterOptions  object
OPP <- parameterOptions()
#Setting new Value for maxSignal
stepSize(OPP)
}
