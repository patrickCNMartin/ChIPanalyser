\name{getHighestFitnessSolutions}
\alias{getHighestFitnessSolutions}

\title{
Get Highest Fitness Solutions
}
\description{
\code{getHighestFitnessSolutions} extract best solution from a ChIPanalyser GA/evolve Run.
}
\usage{
getHighestFitnessSolutions(population,child=2,method="geometric")
}

\arguments{
    \item{population}{
    Population list as output by the \code{evolve} function.
    }
    \item{child}{
    numeric describing the number of solution to be extracted from Population list.
    }
    \item{method}{
    character string describing which scoring method should be used and selected from
    "geometric","ks","MSE","pearson","spearman","kendall",
    "recall","precesion","fscore","MCC","Accuracy" or "AUC".
    }

}
\details{
This function only serves as a way of extracting data from the poppulation list.
Ultimately - it is just a wrapper for some indexing.
}

\value{
Return the index of the top "child" solutions.

}

\author{
Patrick C.N. Martin
}


\examples{



}
