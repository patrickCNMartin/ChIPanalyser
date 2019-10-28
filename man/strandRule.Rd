\name{strandRule}


\alias{strandRule}
\title{
Accessor Method for the \code{strandRule} slot in a
\code{\linkS4class{parameterOptions}} object
}
\description{
Accessor Method for the \code{strandRule} slot in a
\code{\linkS4class{parameterOptions}} object
}
\usage{
strandRule(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\linkS4class{parameterOptions}} object
    }
}
\details{
When computing the PWM Scores and if \code{\link{whichstrand}}
is set to "+-", \code{strandRule} will determine how to handle both strands
( one of three options : "mean", "max", "sum"). If set to "mean",
the average PWM Score of both strand will be computed. If set to "max",
the highest PWM score between each strand will be selected and finally "sum"
will sum both score together.  Default set at "max"
}
\value{
Returns the value assigned to \code{strandRule} slot (one of three options :
"mean", "max", "sum") in a \code{\linkS4class{parameterOptions}} object
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
GPP <- parameterOptions( strandRule="max")
#Accesssing Value for strandRule
strandRule(GPP)
}
