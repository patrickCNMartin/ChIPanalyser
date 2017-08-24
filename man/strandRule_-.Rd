\name{strandRule<-}
\alias{strandRule<-}
\alias{strandRule<--methods}


\title{
Setter method for the \code{strandRule} slot in a
\code{\link{genomicProfileParameters}} object.
}
\description{
Setter method for the \code{strandRule} slot in a
\code{\link{genomicProfileParameters}} object.
}
\usage{
strandRule(object) <- value
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfileParameters}} object
    }
    \item{value}{
    \code{value} is a character string and can be one of the following
    `mean', `max', `sum'. This will only apply if \code{\link{whichstrand}}
    is `+-'. Default set at `max'
    }
}
\details{
When computing the PWM Scores and if \code{\link{whichstrand}} is set
to `+-', \code{strandRule} will determine how to handle both strands
( one of three options : `mean', `max', `sum'). If set to `mean',
the average PWM Score of both strand will be computed. If set to `max',
the highest PWM score between each strand will be selected and finally `sum'
will sum both score together.
Default set at `max'
}
\value{
Returns a \code{\link{genomicProfileParameters}} object with an updated
value for the \code{strandRule} slot
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
GPP <- genomicProfileParameters(PFM=PFM, strandRule="max")
#Setting New Value for strandRule
strandRule(GPP) <- "mean"
}
