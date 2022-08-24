\name{PFMFormat<-}


\alias{PFMFormat<-}

\title{
Setter method for the \code{PFMFormat} slot in a
\code{\link{genomicProfiles}} object
}
\description{
Setter method for the \code{PFMFormat} slot in a
\code{\link{genomicProfiles}} object
}
\usage{
PFMFormat(object) <- value


}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{
    \code{object} is a \code{\link{genomicProfiles}} object
    }
    \item{value}{
    \code{value} is character string of one of the following:
    "raw","transfac","JASPAR" or "sequences". If loading a
    \code{\link{PositionFrequencyMatrix}} from a file, the format
    of the file should specified. Default is raw.
    }
}
\details{
If loading a
\code{\link{PositionFrequencyMatrix}} from a file, the format
of the file should be specified. Default is raw. Please keep in mind that
this argument is used when parsing the
\code{\link{PositionFrequencyMatrix}} file. IF this argument is changed
after building the \code{\link{genomicProfiles}} with
a PositionFrequencyMatrix file, this
will not influence the parsing of the file.
}
\value{
Returns a \code{\link{genomicProfiles}} object with an updated
value for the \code{PFMFormat} slot.
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
# Loading data
data(ChIPanalyserData)
#Loading PFM files
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BEAF-32.pfm")
#Building data objects
#### THIS IS THE PREFFERED METHOD FOR SETTING PFMFormat
GPP <- genomicProfiles(PFM=PFM,PFMFormat="raw")
#Setting New value for PFMFormat
PFMFormat(GPP) <- "raw"
}
