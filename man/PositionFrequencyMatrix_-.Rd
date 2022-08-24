\name{PositionFrequencyMatrix<-}


\alias{PositionFrequencyMatrix<-}

\title{
Setter method for the \code{PFM} slot in a
\code{\link{genomicProfiles}} object
}
\description{
Setter method for the \code{PFM} slot in a
\code{\link{genomicProfiles}} object
}
\usage{
PositionFrequencyMatrix(object)<- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfiles}} object
    }
    \item{value}{
    \code{value} can be of two forms. Either a matrix in the form of a
    Position Frequency Matrix or a path/to/file character string.
    }
}
\details{
The Position Frequency Matrix is one of the fundamental object that needs
to be supplied to  a \code{\link{genomicProfiles}}.
If after building a \code{\link{genomicProfiles}},
only the Position Frequency Matrix needs to be modified then it is
possible to manually update the value of this matrix using the function above.
There are two options for the type of data that may be supplied to the
\code{PFM} slot: a matrix in the form of a Position Frequency Matrix
(matrix with four rows - one for each base pair (ACTG) and a number of
columns equal to the number of sites in the binding site), or it is
possible (also recommended) to provide a path to the file containing the
Position Frequency Matrix. This Position Frequency Matrix file may come
in multiple form such as RAW, Transfac or JASPAR.
WARNING: if a genomicProfiles object has already been created
and only the PFM is supplied/updated ,
then the Positon Weight Matrix will automatically updated as well.
}
\value{
Returns a \code{\link{genomicProfiles}} with an updated \code{PFM}
slot (as described above this will lead to an updated PositionWeightMatrix).
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>}

\examples{
#Loading data
data(ChIPanalyserData)
#Loading PFM files
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BEAF-32.pfm")
#Building genomicProfiles object
GPP<-genomicProfiles()
#Setting PFM
PositionFrequencyMatrix(GPP) <- PFM
}
