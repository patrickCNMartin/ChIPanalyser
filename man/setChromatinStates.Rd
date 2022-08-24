\name{setChromatinStates}
\alias{setChromatinStates}

\title{
setChromatinStates
}
\description{
\code{setChromatinStates} sets chromatin state affinity values to a GRanges object.
}
\usage{
  setChromatinStates(population,chromatinStates)
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
library(ChIPanalyser)
# Input data 
data(ChIPanalyserData)

pop <- 10
params <- c("N","lambda","PWMThreshold", paste0("CS",seq(1:11)))
start_pop <- generateStartingPopulation(pop, params)

cs <- setChromatinStates(start_pop,cs)

}
