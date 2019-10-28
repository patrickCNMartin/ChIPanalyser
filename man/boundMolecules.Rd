\name{boundMolecules}
\alias{boundMolecules}



\title{
Accessor methods for \code{boundMolecules} slot in
\code{\link{parameterOptions}} object.
}
\description{
Extract or Access \code{boundMolecules} slot in
\code{\link{parameterOptions}} object.}
\usage{
boundMolecules(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{parameterOptions}} object.
    }
}
\details{
Defaut value: 1000

When computing occupancy (\code{\link{computeOccupancy}}), a value for the
number of bound Molecules to DNA is needed.
This value can be updated and set in a
\code{\link{parameterOptions}} object.
If the number of molecules is unknown,it is possible to infer this value with
\code{\link{computeOptimal}}.
We strongly encourage to set values when building a
\code{\link{parameterOptions}} object.

}
\value{
Returns \code{boundMolecules} slot in
\code{\link{parameterOptions}} object.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.

}
\author{
Patrick C.N. Martin <pm16057@essex.ac.uk>
}
\examples{
# Building parameterOptions object
OPP <- parameterOptions()
#Checking single value by slot accessor
boundMolecules(OPP)


}
