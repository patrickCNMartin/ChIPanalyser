\name{boundMolecules}
\alias{boundMolecules}
\alias{boundMolecules-methods}


\title{
Accessor methods for \code{boundMolecules} slot in
\code{\link{occupancyProfileParameters}} object.
}
\description{
Extract or Access \code{boundMolecules} slot in
\code{\link{occupancyProfileParameters}} object.}
\usage{
boundMolecules(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{occupancyProfileParameters}} object.
    }
}
\details{
Defaut value: 1000

When computing occupancy (\code{\link{computeOccupancy}}), a value for the
number of bound Molecules to DNA is needed.
This value can be updated and set in a
\code{\link{occupancyProfileParameters}} object.
If the number of molecules is unknown,it is possible to infer this value with
\code{\link{computeOptimal}}.
We strongly encourage to set values when building a
\code{\link{occupancyProfileParameters}} object.

}
\value{
Returns \code{boundMolecules} slot in
\code{\link{occupancyProfileParameters}} object.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.

}
\author{
Patrick C.N. Martin <pm16057@essex.ac.uk>
}
\examples{
# Building occupancyProfileParameters object
OPP <- occupancyProfileParameters()
#Checking single value by slot accessor
boundMolecules(OPP)


}
