\name{chipMean}


\alias{chipMean}


\title{
Accessor method for \code{\link{chipMean}} slot in a
\code{\link{parameterOptions}} object.
}
\description{
Accessor method for \code{\link{chipMean}} slot in a
\code{\link{parameterOptions}} object.

}
\usage{
chipMean(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{parameterOptions}}
    }
}
\details{
Default vlaue : 150
When computing ChIP-seq like profiles (\code{\link{computeChIPProfile}},
the occupancy values given by \code{computeOccupancy} are transformed
into ChIP-seq like profiles.
The average size of a ChIP-seq peak was described by Kaplan
(Kaplan et al. , 2011). It is advised to use the average
width of ChIP peaks from actual ChIP-seq data.
We strongly encourage to set values when building a
\code{\link{parameterOptions}} object.
}

\value{
Returns \code{chipMean} slot from a
\code{\link{parameterOptions}} object.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.

Kaplan T.,Li X.-Y.,Sabo P.J.,Thomas S.,Stamatoyannopoulos J.A., Biggin M.D.,
EisenM.B.
Quantitative models of the mechanisms that control genome-wide patterns of
transcription factor binding during early Drosophila development,
PLoS Genet.,2011, vol. 7 pg. e1001290

}
\author{
Patrick C.N. Martin <pm16057@essex.ac.uk>
}

\examples{
# Building parameterOptions object
OPP <- parameterOptions()
#Accessing chipMean slot in OPP
chipMean(OPP)

}
