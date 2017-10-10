\name{ploidy<-}
\alias{ploidy<-}
\alias{ploidy<--methods}

\title{
Setter Method for the \code{ploidy} slot in an
\code{\link{occupancyProfileParameters}} object
}
\description{
Setter Method for the \code{ploidy} slot in an
\code{\link{occupancyProfileParameters}} object
}

\usage{
ploidy(object)<- value
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{occupancyProfileParameters}} object
    }
    \item{value}{
    \code{value} is a positive integer that describes the ploidy of the
    organism of interest.
    }
}
\details{
Default value for \code{ploidy} is set a 2. It should be mentioned that
ChIPanalyser is based on a model that also considers the ploidy of the
organism however this only considers simple polyploidy (or haploidy). T
he model does not consider hybrids such as wheat.
}
\value{
Returns a \code{\link{occupancyProfileParameters}} object with an updated
value for the \code{ploidy} slot.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>
}
\examples{
# Building occupancyProfileParameters object
OPP <- occupancyProfileParameters()
#Setting new Value for maxSignal
ploidy(OPP) <- 2

}
