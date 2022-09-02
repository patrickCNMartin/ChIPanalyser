\name{PositionWeightMatrix<-}


\alias{PositionWeightMatrix<-}

\title{
Setter Method for the \code{PositionWeightMatrix} slot in a
\code{\link{genomicProfiles}} object
}
\description{
Setter Method for the \code{PositionWeightMatrix} slot in a
\code{\link{genomicProfiles}} object
}
\usage{
PositionWeightMatrix(object) <- value
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfiles}} object
    }
    \item{value}{
    \code{value} is a matrix in the form of a Position Weight Matrix.
    }
}
\details{
If a Position Weight Matrix is readily available, it is possible to directly
assign this matrix to the \code{PWM} slot. However, this is only possible
if a \code{\link{genomicProfiles}} object has already been created.
In that case, we advise to first create a
\code{\link{genomicProfiles}} object. It should be noted
that this Position Weight Matrix will be automatically computed from a
Position Frequency Matrix. If no Position Frequency Matrix are available,
then a Position Weight Matrix can be directly assigned to this slot.
}
\value{
Returns a \code{\link{genomicProfiles}} object with an updated
value for the \code{PWM} slot
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
#Building genomicProfiles object
GPP <- genomicProfiles()
#Setting PWM to PositionWeightMatrix slot
PWM <- matrix(runif(32,-10,20), ncol=8)
rownames(PWM) <- c("A","C","T","G")
PositionWeightMatrix(GPP) <- PWM
}
