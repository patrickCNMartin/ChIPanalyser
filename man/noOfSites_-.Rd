\name{noOfSites<-}
\alias{noOfSites<-}
\alias{noOfSites<--methods}

\title{
Setter Method for the \code{noOfSites} slot in a
\code{\link{genomicProfileParameters}} object.
}
\description{
Setter Method for the \code{noOfSites} slot in a
\code{\link{genomicProfileParameters}} object.

}
\usage{
noOfSites(object) <- value
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfileParameters}} object.
    }
    \item{value}{
    \code{value} is a positive integer that will be assigned to the
    \code{noOfSites} slot.
    }
}
\details{
While computing Position Weight Matricies (PWM) from
Position Frequency Matricies (PFM), it is possible to restrict the number
of sites that will be used to compute the PWM. The default is set at 0.
In this case, all sites will be used to compute the PWM.
}
\value{
Returns a \code{\link{genomicProfileParameters}} object with an updated
value for the \code{noOfSites} slot.
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
GPP <- genomicProfileParameters(PFM=PFM,noOfSites=0)
#Setting new Value for naturalLog
noOfSites(GPP) <- 8

}
