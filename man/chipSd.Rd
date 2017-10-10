\name{chipSd}
\alias{chipSd}
\alias{chipSd-methods}


\title{
Accessor method for \code{chipSd} slot in a
\code{\link{occupancyProfileParameters}} object.
}
\description{
Access or Extract \code{chipSd} slot in a
\code{\link{occupancyProfileParameters}} object.
}
\usage{
chipSd(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{occupancyProfileParameters}}
    }
}
\details{
When computing ChIP-seq like profiles (\code{\link{computeChipProfile}},
the occupancy values given by \code{computeOccupancy} are transformed into
ChIP-seq like profiles.
The average size of a ChIP-seq peak was described by Kaplan
(Kaplan et al. , 2011). The average peak size is subject to
variation. This variation is accounted for with \code{chipSd}.
It is advised to use the standard deviation of ChIP peak width from actual
ChIP-seq data.
We strongly encourage to set values when building a
\code{\link{occupancyProfileParameters}} object.

}
\value{
Returns a \code{\link{occupancyProfileParameters}} object with an
updated value for \code{chipSd}.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.

Kaplan T.,Li X.-Y.,Sabo P.J.,Thomas S.,Stamatoyannopoulos J.A., Biggin M.D.,
Eisen M.B.
Quantitative models of the mechanisms that control genome-wide patterns of
transcription factor binding during early Drosophila development,
PLoS Genet.,2011, vol. 7 pg. e1001290

}
\author{
Patrick C.N. Martin <pm16057@essex.ac.uk>
}

\examples{
# Building occupancyProfileParameters object
OPP <- occupancyProfileParameters()
# Accessing chipSd slot
chipSd(OPP)

}
