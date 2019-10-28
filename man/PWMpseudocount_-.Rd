\name{PWMpseudocount<-}


\alias{PWMpseudocount<-}

\title{
Setter Method for the \code{pseudocount} slot in a
\code{\link{parameterOptions}} object
}
\description{
Setter Method for the \code{pseudocount} slot in a
\code{\link{parameterOptions}} object
}
\usage{
PWMpseudocount(object) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{
    \code{object} is a \code{\link{parameterOptions}} object
    }
    \item{value}{
    \code{value} is a numeric value that will be assigned to the
    \code{pseudocount} slot. Defualt is set at 1
    }
}
\details{
In the context of Position Weight Matricies, the pseudocount is used to
avoid 0 probabilities during the transformation of Position Frequency Matrix
to a Position Probability Matrix and finally to a Postion Weight Matrix.
It is essentially a sample correction that is added in the case of small
sample size. The effect of the base pair to which a pseudocount was
assigned will not influence the model nor will create mathematical issues
such as infinities or zero division.
}
\value{
Returns a \code{parameterOptions} object with an updated value for
the \code{pseudocount} slot.
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
GPP <- parameterOptions( PWMpseudocount=0)
#Setting Value for new PWMpseudocount
PWMpseudocount(GPP) <- 1
}
