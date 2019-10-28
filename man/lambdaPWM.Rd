\name{lambdaPWM}

\alias{lambdaPWM}


\title{
Accessor Method for the \code{lambdaPWM} slot in a
\code{\link{parameterOptions}} object
}
\description{
Accessor Method for the \code{lambdaPWM} slot in a
\code{\link{parameterOptions}} object

}
\usage{
lambdaPWM(object)
}
\arguments{
    \item{object}{
    \code{object} is \code{\link{parameterOptions}} object
}
}
\details{
    The model underlying ChIPanalyser internally infers two paramters: number
    of bound molecules and lambda. Lambda represents a scaling factor for the
    Position weight matrix (PWM). This can be described as how well does a TF
    discriminate between high affinity and very high affinity sites.

}
\value{
Returns the value assigned to the \code{lambdaPWM} slot in a
\code{\link{parameterOptions}} object.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.}
\author{
Patrick C. N. Martin <p.martin@essex.ac.uk>
}

\examples{
# Loading data
data(ChIPanalyserData)

#Building data objects
GPP <- parameterOptions(lambdaPWM=1)
#Setting new Value for lambdaPWM
lambdaPWM(GPP)


}
