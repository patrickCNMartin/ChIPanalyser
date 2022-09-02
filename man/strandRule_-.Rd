\name{strandRule<-}


\alias{strandRule<-}

\title{
Setter method for the \code{strandRule} slot in a
\code{\linkS4class{parameterOptions}} object.
}
\description{
Setter method for the \code{strandRule} slot in a
\code{\linkS4class{parameterOptions}} object.
}
\usage{
strandRule(object) <- value
}

\arguments{
    \item{object}{
    \code{object} is a \code{\linkS4class{parameterOptions}} object
    }
    \item{value}{
    \code{value} is a character string and can be one of the following
    `mean', `max', `sum'. This will only apply if \code{\link{whichstrand}}
    is `+-'. Default set at `max'
    }
}
\details{
When computing the PWM Scores and if \code{\link{whichstrand}} is set
to `+-', \code{strandRule} will determine how to handle both strands
( one of three options : `mean', `max', `sum'). If set to `mean',
the average PWM Score of both strand will be computed. If set to `max',
the highest PWM score between each strand will be selected and finally `sum'
will sum both score together.
Default set at `max'
}
\value{
Returns a \code{\linkS4class{parameterOptions}} object with an updated
value for the \code{strandRule} slot
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
Patrick C.N. Martin and Nicolae Radu Zabe (2020) Dissecting the binding mechanisms of transcription factors to DNA using a statistical thermodynamics framework. CSBJ, 18, 3590-3605.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>
}

\examples{
# Loading data
data(ChIPanalyserData)

#Building data objects
GPP <- parameterOptions(strandRule="max")
#Setting New Value for strandRule
strandRule(GPP) <- "mean"
}
