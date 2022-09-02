\name{removeBackground}


\alias{removeBackground}
\title{
Accessor Method for the \code{removeBackground} slot in a
\code{\link{parameterOptions}} object
}
\description{
Accessor Method for the \code{removeBackground} slot in a
\code{\link{parameterOptions}} object
}
\usage{
removeBackground(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{parameterOptions}} object
    }
}
\details{
A numeric value describing a threshold at which Occupancy signals must be
removed (Default is set at 0). The removal of Occupancy signals will occur
when computing \code{\link{computeOccupancy}}
(see \code{\link{computeOccupancy}} function)
}
\value{
Returns the value assigned to the \code{removeBackground} slot in a
\code{\link{parameterOptions}} object
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
Patrick C.N. Martin and Nicolae Radu Zabe (2020) Dissecting the binding mechanisms of transcription factors to DNA using a statistical thermodynamics framework. CSBJ, 18, 3590-3605.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>

}
\examples{
#Building parameterOptions object
OPP <- parameterOptions()
#Accessing  Value for removeBackground
removeBackground(OPP)

}
