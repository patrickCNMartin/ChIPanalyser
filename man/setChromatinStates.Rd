\name{setchromatinStates}
\alias{setchromatinStates}

\title{
setchromatinStates
}
\description{
\code{setchromatinStates} sets chromatin state affinity values to a GRanges object.
}
\usage{
  setchromatinStates(population,chromatinStates)
}

\arguments{
    \item{population}{
    Population list containing all individuals and associated parameter. Must
    contain chromatin state affinity values. See \code{generateStartingPopulation}.
    }
    \item{chromatinStates}{
    GRanges object containing chromatin state locations.
    }

}
\details{
Chromatin states can be loaded into R as a GRanges object. Each range represents
the extent of a certain chromatin state and the chromatin state type should be
assigned to a meta data column called "name". The affinity values names should be
set accordingly.

}

\value{
Returns a GRange object with affinity scores for each chromatin state range.
Affinity scores are placed in the DNAAffinity meta data column.

}

\author{
Patrick C.N. Martin
}


\examples{



}
