\name{PWMpseudocount}
\alias{PWMpseudocount}
\alias{PWMpseudocount-methods}

\title{
Accessor Method for a \code{PWMpseudocount} slot in a
\code{\link{genomicProfileParameters}}
}
\description{
Accessor Method for a \code{PWMpseudocount} slot in a
\code{\link{genomicProfileParameters}}
}
\usage{
PWMpseudocount(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfileParameters}} object.
    }
}
\details{
In the context of Position Weight Matricies, the pseudocount is used
to avoid 0 probabilities during the transformation of
Position Frequency Matrix to a Position Probability Matrix and finally to a
Postion Weight Matrix. It is essentially a sample correction that is added
in the case of small sample size. The effect of the base pair to which a
pseudocount was assigned will not influence the model nor will create
mathematical issues such as infinities or zero division. Default is set at 1.

}
\value{
Returns the value assigned to a \code{PWMpseudocount} slot in a
\code{\link{genomicProfileParameters}} object
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>

}
\examples{
# Loading data
data(ChIPanalyserData)
#Loading PFM files
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BCDSlx.pfm")
#Building data objects
GPP <- genomicProfileParameters(PFM=PFM, PWMpseudocount=0)
#Accessing slot value
PWMpseudocount(GPP)
}
