\name{naturalLog<-}


\alias{naturalLog<-}

\title{
Setter method for the \code{naturalLog} slot in a
\code{\link{parameterOptions}} object.
}
\description{
Setter method for the \code{naturalLog} slot in a
\code{\link{parameterOptions}} object.
}
\usage{
naturalLog(object)<- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{
    \code{object} is a \code{\link{parameterOptions}} object.
    }
    \item{value}{
    \code{value} is a logical value that will determine if the natural log or
    log2 should be used for the computation of the Position Weight Matrix.
    }
}
\details{
During the computation of a Postion Weight Matrix, the
Position Probability Matrix (derived from a Position Frequency Matrix)
is log transformed. This parameter provides whcih "log transform" will be used.
If TRUE, the Natural Log will bu used (ln). If FALSE, log2 will be used.
We strongly encourage to set values when building a
\code{\link{parameterOptions}} object.
}
\value{
Returns \code{\link{parameterOptions}} object with an updated
value for the \code{naturalLog} slot.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
Patrick C.N. Martin and Nicolae Radu Zabe (2020) Dissecting the binding mechanisms of transcription factors to DNA using a statistical thermodynamics framework. CSBJ, 18, 3590-3605.
}
\author{
Patrick C.N. Martin <p.martin@essex.ac.uk>
}
\examples{
# Loading data
data(ChIPanalyserData)

#Building data objects
OPP <- parameterOptions(naturalLog=TRUE)
#Setting new Value for naturalLog
naturalLog(OPP) <- FALSE

}
