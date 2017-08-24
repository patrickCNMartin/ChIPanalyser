\name{naturalLog}
\alias{naturalLog}
\alias{naturalLog-methods}


\title{
Accessor method the \code{naturalLog} slot in a
\code{\link{genomicProfileParameters}} object.
}
\description{
Accessor method the \code{naturalLog} slot in a
\code{\link{genomicProfileParameters}} object.

}
\usage{
naturalLog(object)
}
\arguments{
    \item{object}{
    \code{object} is \code{\link{genomicProfileParameters}} object.
    }
}
\details{
During the computation of a Postion Weight Matrix, the
Position Probability Matrix (derived from a Position Frequency Matrix)
is log transformed. This parameter provides whcih "log transform" will be used.
If TRUE, the Natural Log will bu used (ln). If FALSE, log2 will be used.
We strongly encourage to set values when building a
\code{\link{genomicProfileParameters}} object.

}
\value{
Returns the value assigned to the \code{naturalLog} slot in a
\code{\link{genomicProfileParameters}} object.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}
\author{
Patrick C.N. Martin <p.martin@essex.ac.uk>
}
\examples{
# Loading data
data(ChIPanalyserData)
#Loading PFM files
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BCDSlx.pfm")
#Building data objects
GPP <- genomicProfileParameters(PFM=PFM,naturalLog=TRUE)
#Setting new Value for naturalLog
naturalLog(GPP)
}
