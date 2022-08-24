\name{getTestingData}
\alias{getTestingData}

\title{
Extract testing data from ChIPscore object
}
\description{
\code{getTestingData} extracts selected regions from ChIPscore object to be used as testing set.
}
\usage{
getTestingData(ChIPscore,loci = 1)
}

\arguments{
    \item{ChIPscore}{
    ChIPscore object as returned by \code{processingChIP}
    }
    \item{loci}{
    numeric describing index of loci to be used as testing data. 
    }
   
}


\value{
Returns ChIPscore object with the selected testing loci.
 
}



\author{
Patrick C.N. Martin <pcnmartin@gmail.com
}


\examples{

library(ChIPanalyser)
data(ChIPanalyserData)
# See GA vignette for usage 
test <- processingChIP(chip,top)
test <- getTestingData(test, 1:2)

}
