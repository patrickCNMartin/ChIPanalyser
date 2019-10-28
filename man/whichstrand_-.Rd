\name{whichstrand<-}


\alias{whichstrand<-}
\title{
Setter method for the \code{whichstrand} slot in a
\code{\linkS4class{parameterOptions}} object
}
\description{
Setter method for the \code{whichstrand} slot in a
\code{\linkS4class{parameterOptions}} object
}
\usage{
whichstrand(object) <- value
}

\arguments{
    \item{object}{
    \code{object} is a \code{\linkS4class{parameterOptions}} object
    }
    \item{value}{
    \code{value} is a character string specifying which strand should be
    used to compute PWM Scores. The three available options are the
    following: "+","-" or "+-". Default is "+-"
    }
}
\details{
PWM Score may be computed on either the positive strand ("+"),
the negative strand ("-") or on both strands ("+-").
}
\value{
Returns a \code{\linkS4class{parameterOptions}} object with an updated
value for the \code{whichstrand} slot
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
whichstrand(GPP) <- "+"
}
