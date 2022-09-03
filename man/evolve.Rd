\name{evolve}
\alias{evolve}

\title{
Running the ChIPanalyser implementation of a Genetic algorithm.
}
\description{
\code{evolve} pushes a starting population to evolve in a genetic algorithm.
}
\usage{
evolve(population,DNASequenceSet,ChIPScore,
       genomicProfiles,parameters=NULL,generations=100,mutationProbability=0.3,
       offsprings=5,chromatinState=NULL,
       method="geometric", lambda=TRUE,
       checkpoint=TRUE,
       filename=NULL, cores=1)
}

\arguments{
    \item{population}{
    numeric value describing the number of individuals in the starting population.
    Alternatively - a starting population list as returned by \code{generateStartingPopulation}.
    NOTE: if numeric - the parameter argument is also required.
    }
    \item{DNASequenceSet}{
    DNAStringSet object containing DNA sequences of interest (Extracted from BSgenome)
    }
    \item{ChIPScore}{
    ChIPScore object  as returned by the \code{processingChIP}function
    }
    \item{genomicProfiles}{
    genomicProfiles object containing minimal information (such as the PWM)
    }
    \item{parameters}{
    vector or list containing each parameter that should be added to the chromosome.
    See \code{generateStartingPopulation}
    }
    \item{generations}{
    numeric describing the number of generation before the Genetic algorithm
    should halt.
    }
    \item{mutationProbability}{
    numeric descrbining the rate of mutations for each surviving individual
    }
    \item{offsprings}{
    numeric descrbining the number of individuals surviving to the next generation
    }
    \item{chromatinState}{
    GRanges object containing chromatin state information. Each state should be
    labled in a meta data column named "name". It is advised to use numeric values
    for each state name.
    }
    \item{method}{
    character string describing the scoring metric that should be used.
    ChIPanalyser offers twelve different metrics: correlation coefficients
     (Pearson, Spearman and Kendall), Mean Squared Error (MSE),
     Kolmogorov–Smirnov Distance, precision, recall, accuracy, F-score,
     Matthew’s correlation coefficient (MCC) and
     Area Under Curve Receiver Operator Characteristic (AUC ROC or just AUC)
    }
    \item{lambda}{
    logical describing if lambda value should be pre-computed. Setting to TRUE
    increases the speed of the algorithm.
    }
    \item{checkpoint}{
    logical describing if population parameters at each generations should be 
    saved. 
    }

    \item{filename}{
    character string that will serve as a prefix to the saved intermediate files.
    }
    \item{cores}{
    numeric describing the number of cores used to run the GA.
    }
}
\details{
ChIPanalyser offers a way of finding optimal solution by using a genetic algorithm.
Instead of running the stadard analysis, TF binding affinities to chromatin states
can be extracted via this more complex method. It should be noted that this
method is better suited for the analysis of chromatin states. While the algorithm
still works with simple DNA Accessibility,
it would potentially take more time for accuracy minor gains.


}

\value{
Returns a named list with three elements.
  \itemize{
  \item{database}{ saves the data frame containing all
  scores for each individual since generation 1}
  \item{population}{ saves the last population with
  chromosome values}
  \item{fitest}{saves the fittest individual
   for a given generation}
  }

}



\author{
Patrick C.N. Martin <pcnmartin@gmail.com
}


\examples{

library(ChIPanalyser)
data(ChIPanalyserData)
# See GA vignette for usage 

}
