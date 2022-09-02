\name{boundMolecules<-}


\alias{boundMolecules<-}


\title{
Setter method for the \code{boundMolecules} slot in a
\code{\link{parameterOptions}} object.
}
\description{
Setter method for the \code{boundMolecules} slot in a
\code{\link{parameterOptions}} object.
}
\usage{
boundMolecules(object)<-value
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{parameterOptions}} object.
    }
    \item{value}{
    \code{value} is a positive integer or vector of positive
    integers describing the number ofmolecules bound to DNA.
    Defaut value is 1000.
    }
}
\details{
Default value: 1000
When computing occupancy (\code{\link{computeOccupancy}}),
a value for the number of bound Molecules to DNA is needed.
This value can be updated and set in a
\code{\link{parameterOptions}} object.
If the number of molecules is unknown,
it is possible to infer this value with \code{\link{computeOptimal}}.
We strongly encourage to set values when building a
\code{\link{parameterOptions}} object.
}
\value{
Returns a \code{\link{parameterOptions}} object with an updated
value for \code{boundMolecules}.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
Patrick C.N. Martin and Nicolae Radu Zabe (2020) Dissecting the binding mechanisms of transcription factors to DNA using a statistical thermodynamics framework. CSBJ, 18, 3590-3605.
}
\author{
Patrick C.N. Martin <pcnmartin@gmail.com>
}

\examples{
# Building parameterOptions object
OPP <- parameterOptions()
# Setting new boundMolecules value in OPP
boundMolecules(OPP) <- 5000
#Checking value in whole object
OPP
#Checking single value by slot accessor
boundMolecules(OPP)

}
