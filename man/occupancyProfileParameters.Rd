\name{occupancyProfileParameters}
\alias{occupancyProfileParameters}

\title{
Occupancy Profile Parameters
}
\description{
\code{occupancyProfileParameters} is a parameter storage object.
Similarly to \code{\link{genomicProfileParameters}},
this object serves the purpose of handling the set of parameters
needed when using certain functions.
}
\usage{

occupancyProfileParameters(ploidy = 2 , boundMolecules = 1000 ,
    backgroundSignal = 0 ,  maxSignal = 1 , chipMean = 150 ,
    chipSd = 150 , chipSmooth = 250 , stepSize = 10 ,
    removeBackground = 0 )
}

\arguments{
    \item{ploidy}{ \code{\link{ploidy}}:
    A numeric Value descibing the ploidy of the organism. Default: 2}
    \item{boundMolecules}{ \code{\link{boundMolecules}}:
    A vector (or single value) containing the number of bound Molecules
    (bound Transcription Factors): Default: 1000}
    \item{backgroundSignal}{ \code{\link{backgroundSignal}}:
    A numeric value descibing the ChIP-seq background Signal
    (average signal from real ChIP seq data). Default: 0}
    \item{maxSignal}{ \code{\link{maxSignal}}:
    A numeric value describing the highest ChIP-seq signal
    (from real ChIP-seq data). Default: 1}
    \item{chipMean}{ \code{\link{chipMean}}:
    A numeric value describing the mean width of a ChIP- seq peak: Default:150}
    \item{chipSd}{ \code{\link{chipSd}}:
    A numeric value describing the standard deviation of ChIP-seq peaks.
    Default: 150 }
    \item{chipSmooth}{ \code{\link{chipSmooth}}:
    A numeric value describing the width of the window used to smooth
    Occupancy profiles into ChIP profiles. Default:250}
    \item{stepSize}{ \code{\link{stepSize}}:
    A numeric value describing the step Size (in base pairs)
    between each ChIP-seq score. Default:10 (Scored every 10 base pairs)}
    \item{removeBackground}{ \code{\link{removeBackground}}:
    A numeric value describing the value at which score should be removed.
    Defualt:0 (If negative scores then remove) }


}
\details{
The purpose of this \code{occupancyProfileParameters} object is to store
paramters used duing certin of the functions available in \code{ChIPanalyser}.
As opposed to \code{genomicProfileParameters}, this object is not required
for a function to run as all paramters have default values ( the object will
be created internally as default values will be used - This also means that
no argument is required for this object to be built). However, it is
strongly advised to set these parameters to increase the goodness of fit
between the predicted Profiles and the real ChIP-seq profiles. As an example,
\code{maxSignal} is the maximum signal of the real ChIP-seq data,
this signal will depend on the data being used thus should be tuned to
opimise the model.

Default values are set at:

ploidy = 2


boundMolecules = 1000


backgroundSignal = 0


maxSignal = 1


chipMean = 150


chipSd = 150


chipSmooth = 250


stepSize = 10


removeBackground = 0


}
\value{
Returns a \code{occupancyProfileParameters} object
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>

}

\seealso{
\code{\link{occupancyProfileParameters-class}}
\code{\link{genomicProfileParameters}}
\code{\link{genomicProfileParameters-class}}
}
\examples{
OPP <- occupancyProfileParameters()
# OR
OPP <- occupancyProfileParameters(ploidy=2,
    boundMolecules=5000,backgroundSignal=0)

}
