\name{genomicProfiles}
\alias{genomicProfiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genomic Profile object
}
\description{
\code{genomicProfiles} is an S4 object serving two purposes: (i) storing internal computed data and (ii) storing paramter options. This object is parsed through the different steps of the pipeline
to facilitate that parsing and changing of paramters.
}
\usage{
genomicProfiles(..., parameterOptions = NULL, genomicProfiles = NULL, ChIPScore = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
  Any of the user available slots in genomicProfiles.
}
  \item{parameterOptions}{
  If some parameters were already previously computed or stored in a parameterOptions,
  parsing this object will use those values instead of the default ones.
}
  \item{genomicProfiles}{
  If some parameters were already previously computed or stored in a genomicProfiles,
  parsing this object will use those values instead of the default ones.
}
  \item{ChIPScore}{
  If some parameters were already previously computed or stored in a ChIPScore,
  parsing this object will use those values instead of the default ones.
}
}
\details{
 The \code{genomicProfiles} object serves the purpose of storing, and parsing paramters and computed data between the different steps of the pipeline. When creating a \code{genomicProfiles}
 object it is possible to use previously computed values by simply parsing the object to the constructor function.
}
\value{
Returns a \code{genomicsProfiles} object with updated slots for all paramters parsed.
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
Patrick C.N. Martin and Nicolae Radu Zabe (2020) Dissecting the binding mechanisms of transcription factors to DNA using a statistical thermodynamics framework. CSBJ, 18, 3590-3605.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>

}

\seealso{
\code{\link{genomicProfiles}}

\code{\link{parameterOptions}}

}
\examples{
PFM <- file.path(system.file("extdata",package="ChIPanalyser"),"BEAF-32.pfm")
genomicProfiles()
genomicProfiles(PFM=PFM,PFMFormat="raw")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
