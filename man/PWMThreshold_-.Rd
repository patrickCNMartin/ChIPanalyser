\name{PWMThreshold<-}


\alias{PWMThreshold<-}
\title{
Setter Method for the \code{PWMThreshold} slot in a
\code{\link{parameterOptions}} object
}
\description{
Setter Method for the \code{PWMThreshold} slot in a
\code{\link{parameterOptions}} object
}
\usage{
PWMThreshold(object) <- value
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{parameterOptions}} object
    }
    \item{value}{
    \code{value} is a numeric value (between 0 and 1) to be assigned to
    the \code{PWMThreshold} slot in
    \code{\link{parameterOptions}} object.
    Default is set at 0.7
    }
}
\details{
The \code{computePWMScore} function requires a so-called PWM Threshold.
This threshold represents the Threshold at which PWM Score should be selected.
The \code{PWMThreshold} is a positive numeric value (between 0 and 1.
If set at 0, all sites will be selected. If set at 0.7 (Default value),
then 70 \% of PWM Score (and by extension binding sites) will be IGNORED.
The top 30 \% will be selected.
}
\value{
Returns  \code{\link{parameterOptions}} objetc with an updated value
for the \code{PWMThreshold} slot
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

#Building data objects
GPP <- parameterOptions(PWMThreshold=0.7)
#Setting Value for new PWMThreshold
PWMThreshold(GPP) <- 0.8
}
