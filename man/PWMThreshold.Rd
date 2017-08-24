\name{PWMThreshold}
\alias{PWMThreshold}
\alias{PWMThreshold-methods}

\title{
Accessor method for the \code{PWMThreshold} slot in  a
\code{\link{genomicProfileParameters}} object
}
\description{
Accessor method for the \code{PWMThreshold} slot in  a
\code{\link{genomicProfileParameters}} object
}
\usage{
PWMThreshold(object)
}

\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfileParameters}} object
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
Returns the value assinged to the \code{PWMThreshold} slot in a
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
GPP <- genomicProfileParameters(PFM=PFM, PWMThreshold=0.7)
#Accessing Value for  PWMThreshold
PWMThreshold(GPP)
}
