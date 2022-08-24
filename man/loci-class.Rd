\name{loci-class}
\Rdversion{1.1}
\docType{class}
\alias{loci-class}


\title{Class \code{"loci"}}
\description{
Setter for Loci of interest parsed to or extracted from the ChIPScore object
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
  \describe{
    \item{.loci<-}{\code{signature(object = "ChIPScore", value = "loci")}: ... }
	 }
}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
Patrick C.N. Martin and Nicolae Radu Zabe (2020) Dissecting the binding mechanisms of transcription factors to DNA using a statistical thermodynamics framework. CSBJ, 18, 3590-3605.
}
\author{
Patrick C. N. Martin <pm16057@essex.ac.uk>

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\linkS4class{ChIPScore}}
}
\examples{
showClass("loci")
}
\keyword{classes}
