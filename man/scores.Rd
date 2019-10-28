\name{scores}

\alias{scores}


\title{
Accessor Method for the \code{scores} slot in a
\code{\linkS4class{ChIPScore}} object
}
\description{
Setter Method for the \code{scores} slot in a
\code{\linkS4class{ChIPScore}} object

}
\usage{
scores(object)
}
\arguments{
    \item{object}{
    \code{object} is \code{\linkS4class{ChIPScore}} object
    }

}
\details{
    When using the \code{\link{processingChIP}}, this functions will return a
    name list of normalised ChIP scores at loci of interest. This functions enalbles
    you to extract those scores from the \code{\linkS4class{ChIPScore} object.}

}
\value{
Returns the value assigned to the \code{scores} slot in a
\code{\linkS4class{ChIPScore}} object.
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

eveLocusChip<-processingChIP(eveLocusChip,eveLocus)
str(scores(eveLocusChip))



}
