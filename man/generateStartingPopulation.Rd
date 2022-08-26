\name{generateStartingPopulation}
\alias{generateStartingPopulation}

\title{
Generate Starting population for ChIPanalyser Genetic algortihm
}
\description{
\code{generateStartingPopulation} generates a starting population with random
traits for each individual
}
\usage{
generateStartingPopulation(population,parameters,names=NULL)
}

\arguments{
    \item{population}{
    numeric value describing the number of individuals in the starting population.
    }
    \item{parameters}{
    vector or list containing each parameter that should be added to the chromosome.
    }
    \item{names}{
    character describing names that should be added to each individual.
    }
}
\details{
\code{generateStartingPopulation} generates a starting poppulation to be used
in the genetic algortihm implemented in ChIPanalyser. There are two main ways a
starting population can be generated:
  \enumerate{
      \item{by name}{
          Using names of each parameter that should be parse to each "chromosome".
          The possible paramters are N, lambda, PWMThreshold, CS
          ( DNAAffinity or DNAAccessibility also works). CS values should also contain
          a numeric value associated to each chromatin state you wish to parse.
          e.g CS1 ... CS14
          This will generate a value by sampling from a set of predefined value
          for each paramters.
      }
      \item{by value range}{
          Using a named list (names for each parameters). Each element of the list
          should contain three numeric values : length of range, min value, max value.
          (Internally - values are parse to \code{runif})
      }
  }
}

\value{
Returns a list of individuals with a random traits

}

\author{
Patrick C.N. Martin
}


\examples{

## by name
param <- c("N","lambda","PWMThreshold","CS1","CS2","CS3")

pop <- generateStartingPopulation(20,param, names = NULL)

# by range 
paramValue <- list(c(10,1,1000),c(10,0,5),c(10,0,0.9),c(10,0,1),c(10,0,1),c(10,0,1))

pop <- generateStartingPopulation(20,paramValue,names= param)

}
