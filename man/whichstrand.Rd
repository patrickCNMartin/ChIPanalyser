\name{whichstrand}


\alias{whichstrand}
\title{
Accessor method for the \code{whichstrand} slot in a
\code{\linkS4class{parameterOptions}} object
}
\description{
Accessor method for the \code{whichstrand} slot in a
\code{\linkS4class{parameterOptions}} object
}
\usage{
whichstrand(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\linkS4class{parameterOptions}} object
    }
}
\details{
PWM Score may be computed on either the positive strand ("+"),
the negative strand ("-") or on both strands ("+-").
}
\value{
Returns on which strand PWM Scores should be computed
( \code{whichstrand} in a \code{\linkS4class{parameterOptions}} object)
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>
}
\examples{
# Loading data
data(ChIPanalyserData)

#Building data objects
GPP <- parameterOptions( whichstrand="+-")
#Setting New Value for whichstrand
whichstrand(GPP)
}
