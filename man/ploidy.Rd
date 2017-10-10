\name{ploidy}
\alias{ploidy}
\alias{ploidy-methods}

\title{
Accessor method for the \code{ploidy} slot in a
\code{\link{occupancyProfileParameters}} object
}
\description{
Accessor method for the \code{ploidy} slot in a
\code{\link{occupancyProfileParameters}} object

}
\usage{
ploidy(object)
}

\arguments{
    \item{object}{
    \code{object} is  a \code{\link{occupancyProfileParameters}} object
    }
}
\details{
Default value for \code{ploidy} is set a 2. It should be mentioned that
ChIPanalyser is based on a model that also considers the ploidy of the
organism of interest however this only considers simple polyploidy
(or haploidy). The model does not consider hybrids such as wheat.

}
\value{
Returns the value assigned to the \code{ploidy} slot in a
\code{\link{occupancyProfileParameters}} object
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
ploidy(OPP)
}
