################################################################################
########################## GA Generic Functions ################################
################################################################################

### Because for some reason i can't make it work otherwise...
### well i could use weird idex systems but this is just cleaner
.initiateGASolutions<-function(solutions){
    strict<-all(sapply(solutions,function(x){all(sapply(x,is.null))}))
    solutions<-lapply(solutions,.initiateSol,strict=strict)
    return(solutions)
}
## Possible soltion generation
## solution is a named list with each of four parameters that should be generated
## if strict is T then this list is an empty named list
## if strict is F then this list is a named list containing min value, max value, number of samples
## for each paramters that needs to be generated


.initiateSol<-function(solutions,strict=TRUE){
    param<-names(solutions)
   
    if(strict){
        solutions<-mapply(.setStrictSoltions,solutions,param,SIMPLIFY = FALSE)
    } else {
        solutions<-mapply(.setfluxRangeSolutions,solutions,param,SIMPLIFY = FALSE)

    }
    return(solutions)
}


## setting Strict Solutions
.setStrictSoltions<-function(solutions,param){
    getName<-param
    #Setting default solutions values
    if(grepl("boundMolecules",getName,ignore.case=T)|
        grepl("Bound molecules",getName, ignore.case=T)|
        grepl("N",getName,ignore.case=T)){
            solutions <- c(1, 10, 20, 50, 100,
                200, 500,1000,1500,2000,2500,3000,4000,5000,10000,20000,50000, 100000,
                200000)
    } else if(grepl("ScalingFactor",getName,ignore.case=T)|
        grepl("Scaling Factor",getName,ignore.case=T)|
        grepl("lambda",getName,ignore.case=T)){
            solutions <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5 ,4 ,4.5, 5)
    } else if(grepl("CS",getName,ignore.case=T)|
        grepl("DNAAccessibility",getName,ignore.case=T)|
        grepl("DNAAffinity",getName,ignore.case=T)){
            solutions <- seq(0,1,by=0.2)
    } else if(grepl("PWMThreshold",getName,ignore.case=T)|
        grepl("PWM Threshold",getName,ignore.case=T)){
            solutions<-seq(0,0.99,by=0.2)
    } else{
        warning("Unused list element in setStrictSoltuons input list")
    }
    return(solutions)
}


## Generating solutions based on range fluctuation
### need to chancge this to seqence instead of runif
## this is to be able to reuse the same soltuon multiple times
.setfluxRangeSolutions<-function(solutions,param){

    getName<-param
    #Setting default solutions values
    if(grepl("boundMolecules",getName,ignore.case=T)|
        grepl("Bound molecules",getName, ignore.case=T)|
        grepl("^N$",getName,ignore.case=T)){
            # getting random values between
            solutions <- round(seq(solutions[1],solutions[2],l = solutions[3]))
    } else if(grepl("ScalingFactor",getName,ignore.case=T)|
        grepl("Scaling Factor",getName,ignore.case=T)|
        grepl("lambda",getName,ignore.case=T)){
            solutions <- round(seq(solutions[1],solutions[2],l = solutions[3]),2)
    } else if(grepl("CS",getName,ignore.case=T)|
        grepl("DNAAccessibility",getName,ignore.case=T)|
        grepl("DNAAffinity",getName,ignore.case=T)){
            solutions <- round(seq(solutions[1],solutions[2],l = solutions[3]),2)
    } else if(grepl("PWMThreshold",getName,ignore.case=T)|
        grepl("PWM Threshold",getName,ignore.case=T)){
            solutions <- round(seq(solutions[1],solutions[2],l = solutions[3]),2)
    } else{
        warning("Unused list element in setStrictSoltuons input list")
    }
    return(solutions)
}

## Creating a population based in random selection of starting solutions
.generateRandomGASolution<-function(population){
    if(is.list(population[[1]])){
    samples<-lapply(population,function(x){res<-lapply(x,sample,size=1)
                                           res$fitness<-NA
                                           res$computed<-FALSE
                                           res$gen<-1
                                           return(res)})
    } else if(is.matrix(population[[1]])) {
    samples<-lapply(population,function(x){res<-as.list(apply(x,2,sample,size=1))
                                           res$fitness<-NA
                                           res$computed<-FALSE
                                           res$gen<-1
                                           return(res)})
    }
    names(samples)<-names(population)
    return(samples)
}


## Geneate starting lambda values

.lambdas <-  function(population,parameters,names=NULL){
    ## first lets build the list squeleton with the appropriate values
    pop<-vector("list",population)

    if(!is.null(names)){names(pop)<-names}

    if(!is.list(parameters)){
        pop<-lapply(pop,function(x){x<-vector("list",length(parameters))
                                names(x)<-parameters
                                return(x)})
    } else{
        pop<-lapply(pop,function(x){x<-parameters})
    }

    # Then lets assign values to each parameter
    pop<-.initiateGASolutions(pop)

    return(pop)
}


## Second option not as light weight but easier to handle
## essentially making less logical comparaisons
.collapseParam <- function(pop){

    if(is.data.frame(pop)){
        popbuf<-pop[,which(colnames(pop)!="fitness" & colnames(pop)!="computed" &
                           colnames(pop)!="gen")]
        param<-apply(popbuf,1,paste0,collapse='',sep='_')


    } else {
        pop<-do.call("rbind",pop)
        popbuf<-pop[,which(colnames(pop)!="fitness" & colnames(pop)!="computed" &
                           colnames(pop)!="gen")]
        param<-apply(popbuf,1,paste0,collapse='',sep='_')

    }
    return(param)
}

## check if solution have already been computed
.IndexOfExistingSolution<-function(currentSolution,database){
    DBbuffer<-.collapseParam(as.data.frame(database))
    PopBuffer<-.collapseParam(currentSolution)

    #idx<-cbind("DB"=which(DBbuffer %in% PopBuffer),"Pop"=which(PopBuffer %in% DBbuffer))
    DBinPop <- match(DBbuffer,PopBuffer)
    PopinDB <- match(PopBuffer, DBbuffer)

    idx<-list("DB"=DBinPop,"Pop"=PopinDB)
    return(idx)
}

## generation tracking
.fittyMacFitFace<-function(population,method,fit,generation){
  ## extracting fitness and filtering fitness
  Fit<-sapply(population,function(x){x$fitness})
  indexFit<-which(!is.na(Fit))

  ## population template
  bufferPopulation<-population[indexFit]

  ## extracting top solutions
  wellFit <- getHighestFitnessSolutions(bufferPopulation,child=1,method=method)
  FitPopulation<-bufferPopulation[wellFit]

  fit[[generation]]<-FitPopulation[[1]]$fitness
  return(fit)
}


## building database of pre computed solutions
.buildDatabase<-function(population,database,new=FALSE, generation=1){
    if(new==FALSE){

       #if(length(population)!=10)browser()
        idx<-.IndexOfExistingSolution(population,database)


        # Population that needs to be computed
        # Not present in database
        PopinDB<-which(is.na(idx[[2]]))

        ToBeComputedPop <-population[PopinDB]

            # Population that is present in database
            # reassign score
            #DbinPop<- idx[[1]][which(!is.na(idx[[1]]))]
            popbuffer <- idx[[2]][which(!is.na(idx[[2]]))]
            ComputedPop <- population[which(!is.na(idx[[2]]))]




            for(i in seq_along(ComputedPop)){
                ComputedPop[[i]][["fitness"]]<-unlist(database[popbuffer[i],"fitness"])
                ComputedPop[[i]][["computed"]]<-TRUE
                ComputedPop[[i]][["gen"]]<-generation

            }

        return(list("DataBase"=as.data.frame(database),"ComputedPop"=ComputedPop,
                    "ToBeComputedPop"=ToBeComputedPop))
    } else {
        #population <- lapply(population,unlist)
        population<-as.data.frame(do.call("rbind", population))
       
        population$gen<-rep(generation,nrow(population))
        database<-rbind(database,population)
       
        return(database)
    }


}

## Alright lets get to mutation stuff shall we.

.mutate <- function(population,parameters){
    ## generate tempalte of values
    pop<-vector("list",1)
    if(!is.list(parameters)){
        pop<-lapply(pop,function(x){x<-vector("list",length(parameters))
                                names(x)<-parameters
                                return(x)})
    } else{
        pop<-lapply(pop,function(x){x<-parameters})
    }
    solutions<-.initiateGASolutions(pop)
    
    ## mutating individuals
    # where and what to mutate
    mutation<-sample(seq_along(solutions[[1]]),size=1)

    # mutating individuals
    population[[mutation]]<-sample(solutions[[1]][[mutation]][which(solutions[[1]][[mutation]]!=population[[mutation]])],size=1)

    # setting fitness and computed status to default values
    population$fitness<-NA
    population$computed<-FALSE

    # returning New population even if it is called old population
    return(population)
}


### let us do some cross over shit
## radu solutions seems to be very much
## based on random mixing but what iff we take highest
## nope that wont work
## you don't know which parameters is the strongets
## could actually try and do this though (it would super heavy computationanly)
## thats not how you spell that word dipshit

.crossover<-function(ga1, ga2, extent="random"){

    ## selecting the extent of the crossover
    if(extent=="random"){
        # generating random crossover extent
        extent<-sample(seq(from=1,to=(length(ga1)-1)),size=1)

        #selection of cross over regions i guess
        crossoverSol <- sample(seq_along(ga1),size=extent)


    } else {
        # selecting cross over regions with provided extent
        if(extent<(1/length(ga1))){stop("Extent argument out of bounds - Limited by Number of parameters")}
        crossoverSol <- sample(seq_along(ga1),size=floor(length(ga1)*extent))

    }

    #crossing over some stuff on both solutions
    ## using loops because you suck at R
    for(ga in crossoverSol){
        ga1[[ga]]<-ga2[[ga]]
    }

    # New solution with cross over
    ga1$fitness<-NA
    ga1$computed<-FALSE

    return(ga1)
}


## Alright final step

### Generating a new population

.generateNewPopulation <- function(population,parameters,mutationProbability=0.3,offsprings=1,method="geometric",gen=2){
    ## Generate template for new population
    NewPop<-vector("list",length(population))

    ## extracting fitness and filtering fitness
    Fit<-sapply(population,function(x){x$fitness})
    indexFit<-which(!is.na(Fit) & !duplicated(population))

    ## population template
    bufferPopulation<-population[indexFit]

    ## setting number of offsprings
    if(offsprings > floor(length(bufferPopulation)/2) | offsprings <= 0){
        offsprings <- floor(length(bufferPopulation)/2)
    }

    ## extracting top solutions
    wellFit <- getHighestFitnessSolutions(bufferPopulation,child=offsprings,method=method)
    FitPopulation<-bufferPopulation[wellFit]

    ## keeping highest ranking individuals
    for(i in seq_along(wellFit)){
        NewPop[[i]]<-FitPopulation[[i]]
    }

    ## selecting the other solutions that will be use
    indexGAReuse <- seq(length(wellFit)+1,length(NewPop))

    reuseGA <- sample(seq_along(FitPopulation),size=length(indexGAReuse),replace=TRUE)

    # this is just assinging new values for each indiv
    # random sampling of the fit population
    # I feel like there could be a local minimum problem here
    # it's Either mutate OR crossover
    # Should we include possibility of both
    for(j in seq_along(indexGAReuse)){
        localPopulation<-FitPopulation[[reuseGA[j]]]
        ## Taking reused population and mutating it
        if(runif(1)<=mutationProbability){
            
            localPopulation<-.mutate(localPopulation,parameters)
        } else {
             
            localPopulation<-.crossover(localPopulation,NewPop[[sample(seq_along(wellFit),size=1)]])
            #localPopulation<-.crossover(localPopulation,NewPop[[sample(seq_len(j),size=1)]])
        }
        NewPop[[indexGAReuse[j]]]<-localPopulation
    }

    ## return New population
    ## It will probably be easy not to check if in database now.
    # that can just be part of the analysis function
    ## not sure how much of a good idea this is
    ## this might not be released for ChIPanalyser so we can do that later

    return(NewPop)

}







## Single individual Run
## computing the fitness of an individual

.individualSurvival <- function(indiv,DNAAffinity,
    genomicProfiles,DNASequenceSet,ChIPScore,fitness="geometric",lambdaDB=NULL){
    ##### you will definitly need some checks here
    if(!indiv$computed ){

        # Setting lambda from GA population
        lambdaPWM(genomicProfiles)<-unlist(indiv[
            grepl("ScalingFactor",names(indiv),ignore.case=T)|
            grepl("Scaling Factor",names(indiv),ignore.case=T)|
            grepl("lambda",names(indiv),ignore.case=T)])

        # Setting PWMThreshold from GA population
        #PWMThreshold(genomicProfiles)<-unlist(indiv[
            #grepl("PWMThreshold",names(indiv),ignore.case=T)|
            #grepl("PWM Threshhold",names(indiv),ignore.case=T)])
        # Setting Number of bound molecules from GA population

        boundMolecules(genomicProfiles)<- unlist(indiv[
            grepl("boundMolecules",names(indiv),ignore.case=T)|
            grepl("Bound molecules",names(indiv),ignore.case=T)|
            grepl("^N$",names(indiv),ignore.case=T)])

        if(is.null(lambdaDB)){
        ### compitng prediction from GA population
        ## genome wide scoring
        genomeWide <- computeGenomeWideScores(genomicProfiles,DNASequenceSet
            ,DNAAffinity,verbose=FALSE)
        } else{

        lambda<-as.vector(as.matrix(lambdaDB[which(lambdaDB[,1]==lambdaPWM(genomicProfiles)),]))
        genomeWide <- genomicProfiles
        .maxPWMScore(genomeWide)<-lambda[2]
        .minPWMScore(genomeWide)<-lambda[3]
        .averageExpPWMScore(genomeWide)<-lambda[4]
        .DNASequenceLength(genomeWide) <- lambda[5]
        }


        ## compute PWMscores

        ## check naming of pwmscore its a bit shit also you need to find a more elegant way to deal with it
        pwmScores <- computePWMScore(genomeWide,DNASequenceSet,ChIPScore,DNAAffinity, verbose=FALSE)
        if(is.null(pwmScores)){
            return(indiv)
        }


        ## compute occupancy
        occup <- computeOccupancy(pwmScores,verbose=FALSE)

        ## computing ChIP profiles
        chip <- computeChIPProfile(occup,ChIPScore,verbose=FALSE)

        ## setting the right method
        if(length(grep(fitness,c("fscore","AUC","recall","precision","accuracy","MCC"),ignore.case=TRUE))>0){
            internalMethod<-"fscore"
        } else {
            internalMethod <- fitness
        }

        ## Goodness of fit

        GoF <- profileAccuracyEstimate(chip,ChIPScore,method=internalMethod)

        ## extract top scores

        #topScore <- GoF[[1]][[1]][[1]][paste0(fitness,"Mean")]
        topScore <- profiles(GoF)[[1]][[1]][paste0(fitness,"Mean")]
        ## updating indiv
        indiv$fitness <- topScore
        indiv$computed <- TRUE


    }
    # return idividual with updated fitness
    return(indiv)

}




### extracting values from databse
.populationFitness<- function(database,parameters){
    ## extract values
    database<-as.data.frame(database)
    fitness <- unlist(database$fitness)
    param <- .collapseParam(database)

    ## splitting parameters
    param<-strsplit(param,"_")

    ## building matrix
    optimal<-.lambdas(1,parameters)[[1]]
    ## this will need to be more robust to ensure that you have name flexibility
    mat<-matrix(0,ncol=length(optimal[["N"]]),nrow=length(optimal[["lambda"]]))
    colnames(mat)<-optimal[["N"]]
    rownames(mat)<-optimal[["lambda"]]
    ## replacing values in matrix
    for(i in seq_along(fitness)){
       if(!is.na(fitness[i])){
            rows<-which(rownames(mat)==param[[i]][2])
            cols<-which(colnames(mat)==param[[i]][1])
            mat[rows,cols]<-fitness[i]
       } else{
            rows<-which(rownames(mat)==param[[i]][2])
            cols<-which(colnames(mat)==param[[i]][1])
            mat[rows,cols]<-0
       }
    }

    return(mat)
}


## chunk data for training a validation 

getTrainingData <- function(ChIPscore,loci = 1){
    train <- ChIPscore
    .scores(train) <- scores(train)[loci]
    .loci(train) <- loci(train)[loci]
   ## quick check to make sure it's not empty
    if(length(scores(train))<1){
        stop("Empty training set!")
    }
    return(train)
}

getTestingData <- function(ChIPscore,loci = 1){
    validation <- ChIPscore
    .scores(validation) <- scores(validation)[loci]
    .loci(validation) <- loci(validation)[loci]
    
    if(length(scores(validation))<1){
        stop("Empty validation set!")
    }
    return(validation)
}

splitData <- function(ChIPscore, dist = c(80,20), as.proportion = TRUE){
    locs <- length(scores(ChIPscore))
    if(locs <= 1){
        stop("Cannot split less than 2 loci!")
    }
    if(as.proportion){
        train <- floor(locs*(dist[1]/100))
        if(train < 1)stop("Training set proportion too small!")
        test <- floor(locs*(dist[1]/100)) + floor(locs*(dist[2]/100))
        if(test < 1)stop("Validation set proportion too small!")

        trainSet <- getTrainingData(ChIPscore, loci = seq(1,train))
        testSet <- getTestingData(ChIPscore,  loci = seq(train +1 , test))
    }else {
        if(length(dist) != 4){
            stop("Please Provide start/stop loci index for training and testing")
        }
        trainSet <- getTrainingData(ChIPscore, loci = seq(dist[1],dist[2]))
        testSet <- getTestingData(ChIPscore,  loci = seq(dist[3], dist[4]))
    }
   

    return(list("trainingSet" = trainSet, "testingSet" = testSet))

}