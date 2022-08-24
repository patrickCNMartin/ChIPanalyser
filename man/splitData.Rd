\name{splitData}
\alias{splitData}

\title{
Get Training and Testing data from ChIPscore objects 
}
\description{
\code{splitData} splits processed ChIP data into training and testing sets. 
}
\usage{
splitData(ChIPscore,dist = c(80,20), as.proportion = TRUE)
}

\arguments{
    \item{ChIPscore}{
    ChIPscore object as returned by \code{processingChIP}
    }
    \item{dist}{
    If \code{as.proportion} is to TRUE, split the data into 
    desired proportions. Default sets 80\% training and 
    20\% testing.
    If \code{as.proportion} is to FALSE,a vector of 4 numeric 
    values describing start and end of training and testing 
    respectively. 
    }
    \item{as.proportion}{
    Logical describing if values provided to \code{dist} should
    be treated as \% of training and testing or if \code{dist}
    should be  considered as start and end of loci selected for 
    training and testing respectively.
    }
   
}


\value{
Returns a named list of ChIPScore objects

* trainingSet = ChIPscore containing training set
* testingSet = ChIPscore containing testing set.  
}



\author{
Patrick C.N. Martin <pcnmartin@gmail.com
}


\examples{

library(ChIPanalyser)
data(ChIPanalyserData)
# See GA vignette for usage 
test <- processingChIP(chip,top)
usingDist <- splitData(test, dist = c(50,50),as.proportion = TRUE )
usingIndex <- splitData(test, dist = c(1,2,3,4),as.proportion = FALSE )

}
