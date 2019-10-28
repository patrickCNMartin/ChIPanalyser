\name{noiseFilter<-}

\alias{noiseFilter<-}


\title{
Setter Method for the \code{noiseFilter} slot in a
\code{\link{parameterOptions}} object
}
\description{
Setter Method for the \code{noiseFilter} slot in a
\code{\link{parameterOptions}} object

}
\usage{
noiseFilter(object) <- value
}
\arguments{
    \item{object}{
    \code{object} is \code{\link{parameterOptions}} object
    }
    \item{value}{\code{value} is the value to be assigned to the
    \code{noiseFilter} slot (zero - mean - median -sigmoid)}
}
\details{
Noise filtering method that
should be used on ChIP-seq data. Four methods are available: Zero, Mean,
Median and Sigmoid. Zero removes all ChIP-seq scores bellow zero, mean under
the mean score, median under median score and sigmoid assignes a weight
to each score based on a logistic regression curve. Mid point is set at 95
95 quantile of ChIP-seq scores. Below midpoint will receive a score between
0 and 1 , everything above will receive a score between 1 and 2

}
\value{
Returns the value assigned to the \code{noiseFilter} slot in a
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
GPP <- parameterOptions(noiseFilter="sigmoid")
#Setting new Value for noiseFilter
noiseFilter(GPP) <-"zero"


}
