\name{noOfSites}

\alias{noOfSites}


\title{
Accessor Method for the \code{noOfSites} slot in a
\code{\link{parameterOptions}} object
}
\description{
Accessor Method for the \code{noOfSites} slot in a
\code{\link{parameterOptions}} object

}
\usage{
noOfSites(object)
}
\arguments{
    \item{object}{
    \code{object} is \code{\link{parameterOptions}} object
}
}
\details{
While computing Position Weight Matricies (PWM) from
Position Frequency Matricies (PFM), it is possible to restrict the
number of sites that will be used to compute the PWM.
The default is set at "all". In this case, all sites will be used to
compute the PWM.

}
\value{
Returns the value assigned to the \code{noOfSites} slot in a
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
GPP <- parameterOptions(noOfSites="all")
#Setting new Value for naturalLog
noOfSites(GPP)


}
