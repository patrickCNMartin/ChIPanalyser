\name{chipMean<-}
\alias{chipMean<-}
\alias{chipMean<--methods}

\title{
Access methods for \code{chipMean} slot in
\code{\link{occupancyProfileParameters}} object.
}
\description{
Access methods for \code{chipMean} slot in
\code{\link{occupancyProfileParameters}} object.
}
\usage{
chipMean(object)<-value
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{occupancyProfileParameters}} object.}
    \item{value}{
    \code{value} is a positive numeric value that will be assigned to the
    \code{chipMean} slot.
    \code{chipMean} describes the average size of a
    ChIP-seq peak in base pairs.
    }
}
\details{
Default vlaue : 150
When computing ChIP-seq like profiles (\code{\link{computeChipProfile}},
the occupancy values given by \code{computeOccupancy} are
transformed into ChIP-seq like profiles.
The average size of a ChIP-seq peak was described by Kaplan
(Kaplan et al. , 2011). It is advised to use the average
width of ChIP peaks from actual ChIP-seq data.
We strongly encourage to set values when building a
\code{\link{occupancyProfileParameters}} object.

}
\value{
Returns a \code{\link{occupancyProfileParameters}} object with an updated
value for \code{chipMean} slot.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.

Kaplan T.,Li X.-Y.,Sabo P.J.,Thomas S.,Stamatoyannopoulos J.A.,
Biggin M.D.,EisenM.B.
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
# Setting new value for slot
chipMean(OPP) <- 250
}
