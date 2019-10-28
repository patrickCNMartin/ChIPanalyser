\name{ChIPanalyserData}
\alias{ChIPanalyserData}
\alias{Access}

\alias{eveLocus}
\alias{eveLocusChip}
\alias{geneRef}
\docType{data}
\title{
ChIPanalyserData
}
\description{
\code{ChIPanalyserData} is derived from real biological data.
The source organism is \emph{Drosophila melanogaster}.
The data can be described as genomic data as it contains DNA sequences,
loci, genetic information, DNA accessibility data and ChIP-seq data.
}
\usage{
    data(ChIPanalyserData)
}
\format{
\enumerate{

    \item{\code{Access}}{is \code{\link[GenomicRanges]{GRanges}}
    containing DNA Accesibility data for the sequences described above.}
    \item{\code{eveLocus}}{is \code{\link[GenomicRanges]{GRanges}}
    containing a locus of interest.
    In this case \emph{eve strip Locus} on
    chromosome 2R in \emph{Drosophila melanogaster}}
    \item{\code{eveLocusChip}}{is a GRanges containing ChIP score of the eve
    strip locus in \emph{Drosophila melanogaster}. }
    \item{\code{geneRef}}{is a \code{\link[GenomicRanges]{GRanges}}
    containing genetic information (exon, intron, 3'UTR, 5'UTR)
    for the sequence described above.}
    }

}
\value{
Returns a set of Rdata objects as described above.
}
\source{
Transcription Factor PFM:
Berkeley Drosophila Transcription Network Project (bdtnp.lbl.gov)

}
\references{
Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.
}

\examples{
data(ChIPanalyserData)
}
