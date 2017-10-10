\name{NoAccess}
\alias{NoAccess}
\alias{NoAccess-methods}

\title{
Accessor Method for the \code{NoAccess} slot in a
\code{\link{genomicProfileParameters}} object.
}
\description{
Accessor Method for the \code{NoAccess} slot in a
\code{\link{genomicProfileParameters}} object.
}
\usage{
NoAccess(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfileParameters}} object.
    }
}
\details{
During certain computations, it is possible that the Loci of interest
do no show any overlap with accesible DNA. If this were to be the case,
a warning message will appear in the console but these inaccessible Loci
will be stored in this slot. It is also for these reasons that it is
imperative for Loci of interest to be named
(in this case, a named \code{\link[GenomicRanges]{GRanges}}).
}
\value{
Returns a character string with loci containing no accesible DNA.
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
GPP <- genomicProfileParameters(PFM=PFM)

# Loci with no acces - a warning message will be issued
#if loci do no contain accesible DNA
# Otherwise this slot will remain empty

NoAccess(GPP)
}
