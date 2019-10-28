\name{lociWidth}

\alias{lociWidth}


\title{
Accessor Method for the \code{lociWidth} slot in a
\code{\link{parameterOptions}} object
}
\description{
Setter Method for the \code{lociWidth} slot in a
\code{\link{parameterOptions}} object

}
\usage{
lociWidth(object)
}
\arguments{
    \item{object}{
    \code{object} is \code{\link{parameterOptions}} object
    }

}
\details{
    When using the \code{\link{processingChIP}} function, the provided ChIP
    scores will be split into bins of a given size. lociWidth determines the Size
    of that bin. Default is set at 20 000 bp.
    This means that the ChIP profiles provided will be split into bins of 20 000 bp
    over the entire profile provided if no loci of interest is provided.

}
\value{
Returns the value assigned to the \code{lociWidth} slot in a
\code{\link{parameterOptions}} object.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.}
\author{
Patrick C. N. Martin <p.martin@essex.ac.uk>
}

\examples{
# Loading data
data(ChIPanalyserData)

#Building data objects
GPP <- parameterOptions(lociWidth=20000)
#Accessing new Value for lociWidth
lociWidth(GPP)


}
