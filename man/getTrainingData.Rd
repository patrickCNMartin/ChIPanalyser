\name{getTrainingData}
\alias{getTrainingData}

\title{
Extract training data from ChIPscore object
}
\description{
\code{getTrainingData} extracts selected regions from ChIPScore object to be used as training set.
}
\usage{
getTrainingData(ChIPscore,loci = 1)
}

\arguments{
    \item{ChIPscore}{
    ChIPscore object as returned by \code{processingChIP}
    }
    \item{loci}{
    numeric describing index of loci to be used as training data. 
    }
   
}


\value{
Returns ChIPscore object with the selected training loci.
 
}



\author{
Patrick C.N. Martin <pcnmartin@gmail.com
}


\examples{

library(ChIPanalyser)
data(ChIPanalyserData)
# See GA vignette for usage 
test <- processingChIP(chip,top)
test <- getTrainingData(test, 1:2)

}
