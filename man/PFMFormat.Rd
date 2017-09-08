\name{PFMFormat}
\alias{PFMFormat}
\alias{PFMFormat-methods}

\title{
Accesor method for the \code{PFMFormat} slot in a
\code{\link{genomicProfileParameters}} object
}
\description{
Accesor method for the \code{PFMFormat} slot in a
\code{\link{genomicProfileParameters}} object
}
\usage{
PFMFormat(object)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfileParameters}} object
    }

}
\details{
If loading a
\code{\link{PositionFrequencyMatrix}} from a file, the format
of the file should be specified. Default is raw. Please keep in mind that
this argument is used when parsing the
\code{\link{PositionFrequencyMatrix}} file. IF this argument is changed
after building the \code{\link{genomicProfileParameters}} with
a PositionFrequencyMatrix file, this
will not influence the parsing of the file.
\code{PFMFormat} can be one of the following:
"raw","transfac","JASPAR" or "sequences"
}
\value{
Returns the value assigned to the \code{PFMFormat} slot a
\code{\link{genomicProfileParameters}}
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
#### THIS IS THE PREFFERED METHOD FOR SETTING PFMFormat
GPP <- genomicProfileParameters(PFM=PFM,PFMFormat="raw")
#Setting New value for PFMFormat
PFMFormat(GPP)
}
