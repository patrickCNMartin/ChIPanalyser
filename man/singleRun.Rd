\name{singleRun}
\alias{singleRun}

\title{
singleRun
}
\description{
\code{singleRun} runs ChIPanalyser after optimal paramters have been found by the
evolve function.
}
\usage{
  singleRun(indiv,DNAAffinity,
            genomicProfiles,DNASequenceSet,
            ChIPScore,fitness="geometric")
}

\arguments{
    \item{indiv}{
    Population list containing the top scoring individual. Note that this should be
    a list of length 1 containing another list with all parameter values.
    }
    \item{DNAAffinity}{
    GRanges object as outputed by the \code{setchromatinStates}.
    }
    \item{genomicProfiles}{
    genomicProfiles object containing PWM scores and other desired metrics. Note
    that PWMThreshold, lambda and N will be overwritten using values from indiv.
    }
    \item{DNASequenceSet}{
    DNA string set object containing DNA sequence of interest.
    }
    \item{ChIPScore}{
    ChIPScore object as outputed by the \code{processingChIP} function.
    }
    \item{fitness}{
    character string describing which metric should be used to assess fitness and
    should be one of the following:"geometric","ks","MSE","pearson","spearman","kendall",
    "recall","precesion","fscore","MCC","Accuracy" or "AUC".
    }

}
\details{
Once the genetic algorithm has been optimised, the top individual may be run
on its own to get predicted ChIP profiles. The use of this function requires a
few extract steps in order to predict ChIP profiles.

First, the index of the top individual should be extracted (see \code{getHighestFitnessSolutions}).
Second, using this index, subset top individual from GA population. Note this
should be done using "[]" single bracket notation as, a list of length 1 containing
another list with all parameter values is required for the next steps. Yes, this
is might seem annoying but the functions were design for list structures...
Third, setchromatinStates using the top individual list. This will add chromatin
affinity values to your chromatinState GRanges. Use this new chromatinState object
as your new chromatinState object.
Fourth, parse your indiv list object to \code{singleRun}.


}

\value{
Return a list with three elements.
First element contains a genomicProfiles object with occupancy scores.
Second element contains a genomicProfiles objecy with ChIP profile scores.
Third element contains a goodness of fit metrics.

}

\author{
Patrick C.N. Martin
}


\examples{



}
