##################################################
########## Analaysis function ####################
##################################################



#### Running GA For n generations


evolve <- function(population,DNASequenceSet,ChIPScore,
                   genomicProfiles,parameters=NULL,generations=100,mutationProbability=0.3,
                   offsprings=5,chromatinState=NULL,
                   method="geometric", lambda=TRUE,checkpoint=TRUE,
                   filename=NULL, cores=1){

    ## setting up population
    if(is.numeric(population) & !is.null(parameters)){
      population <- generateStartingPopulation(population,parameters)
    }


      ## setting up chromatinState from starting population
    if(!is.null(chromatinState) & length(grep("cs",names(population[[1]]),ignore.case=T))>1){
        CS <- setChromatinStates(population,chromatinState)
    } else if(!is.null(chromatinState) & length(grep("cs",names(population[[1]]),ignore.case=T))==1){
        CS<-list(chromatinState)
    } else{
        CS<-vector("list", length(population))
    }

    ## building starting database
    database <- as.data.frame(matrix(0,ncol=length(population[[1]])))

    colnames(database)<-names(population[[1]])


    ## setting lambda data based
    if(lambda==TRUE){
        message("Generating Lambda DataBase")
        lambdas <- .lambdas(1,c("lambda"))[[1]]
        lambdaPWM(genomicProfiles)<-unlist(lambdas)
        lambdaDB <- computeGenomeWideScores(genomicProfiles,DNASequenceSet,CS[[1]], cores=cores, verbose=FALSE)
        lambdaDB <-cbind(lambdaPWM(lambdaDB),
            rep(maxPWMScore(lambdaDB),length(lambdas)),
            rep(minPWMScore(lambdaDB),length(lambdas)),
            averageExpPWMScore(lambdaDB),
            rep(DNASequenceLength(lambdaDB),length(lambdas)))
    }

    ## checkpoint saves
    if(checkpoint & lambda){
        if(!is.null(filename)){
            save(lambdaDB,file=paste0(filename,"_",method,"_LambdaDBTraining.Rda"))
        } else {
            save(lambdaDB,file=paste0(method,"_LambdaDBTraining.Rda"))

         }
    }
    ### vec for fittyMacFitFace
    GenFit<-vector("list", generations)

    ## time for some evolution
    for(i in seq_len(generations)){

        message(paste0("Generation: ",i))
        # checking if solution exists or building databe on first generation

        buffer <- .buildDatabase(population,database,generation=i)
        database <- buffer[[1]]
        # removing duplicated enteries when necessary.. we are not going to bother with the computed ones
        #if(length(buffer[[2]])!=0){
          #ComputedPopulation <- buffer[[2]][!duplicated(.collapseParam(buffer[[2]]))]
        #} else{
          ComputedPopulation <- buffer[[2]]
        #}

        ToBeComputedPop <- buffer[[3]]
        




        if(length(ToBeComputedPop)!=0){
          ## re set Chromatin States
          if(!is.null(chromatinState) & length(grep("cs",names(ToBeComputedPop[[1]]),ignore.case=T))>1){
              CS <- setChromatinStates(ToBeComputedPop,chromatinState)
          } else if(!is.null(chromatinState) & length(grep("cs",names(ToBeComputedPop[[1]]),ignore.case=T))==1){
              CS<-list(chromatinState)
          } else{
              CS<-vector("list", length(ToBeComputedPop))
          }
        # Compute fitness of idividual

        if(lambda==TRUE){

            ToBeComputedPop <- parallel::mcmapply(.individualSurvival,ToBeComputedPop,CS,
            MoreArgs=list(genomicProfiles,DNASequenceSet,
            ChIPScore=ChIPScore,fitness=method,lambdaDB=lambdaDB),SIMPLIFY=FALSE,mc.cores=cores)
        } else {
            ToBeComputedPop <- parallel::mcmapply(.individualSurvival,ToBeComputedPop,CS,
            MoreArgs=list(genomicProfiles,DNASequenceSet,
            ChIPScore=ChIPScore,fitness=method,lambdaDB=NULL),SIMPLIFY=FALSE,mc.cores=cores)
        }


        # update database

        database <- .buildDatabase(ToBeComputedPop,database,new=TRUE,generation=i)
        database <- database[!duplicated(.collapseParam(database)),]
        }
        # Rebuilding population with both "pre/post" computed

        population<-c(ComputedPopulation,ToBeComputedPop)

        # Best performance
        GenFit<-.fittyMacFitFace(population,method,GenFit,i)




        message(paste0("Generation Fitness: ",GenFit[[i]]))
        ## checkpoint saving this is messy af
        # shame shame shame Patrick
        if(checkpoint){
            if(!is.null(filename)){
                fileCheckPop<-paste0(filename,"_",method,"_population.csv")
             }else{
                fileCheckPop<-paste0(filename,"_",method,"_population.csv")
             }
            if(i>1){
                out <- .revert_to_class(do.call("rbind",population))
                write.table(out,
                    file = fileCheckPop,
                    sep = ",",
                    col.names= TRUE, 
                    row.names = FALSE,
                    append = TRUE )
            } else {
                out <- .revert_to_class(do.call("rbind",population))
                write.table(out,
                    file = fileCheckPop,
                    sep = ",",
                    col.names= TRUE, 
                    row.names = FALSE,
                    append = FALSE)
            }

        }

        # create new population
        if(i != generations){
        ## check for odd stuff

        population <- .generateNewPopulation(population,parameters,
            mutationProbability=mutationProbability,
            offsprings=offsprings, method=method,gen=i)



        }
    }

    ## When generation are over return databse and last population
    return(list("database"=as.data.frame(apply(database,2,unlist)),
                "population"=population,
                "fitest"=GenFit))

}

getHighestFitnessSolutions<-function(population,child=2,method="geometric"){

    ## checking if you have the right method
    if(length(grep(method,c("geometric","ks","MSE","pearson","spearman","kendall",
    "recall","precesion","fscore","MCC","Accuracy","AUC"),ignore.case=T))<1){
        stop("Selected method is not part of the following:
        geometric,ks,MSE,pearson,spearman,kendall,
        recall,precesion,fscore,MCC,Accuracy,AUC")
    }


    ## We are just going to consider a interger value for children. we can add that option somewhere else
    ## in the create new pop function
        
        UFit<-unlist(sapply(population,"[[","fitness"))

        if(grepl("ks",method,ignore.case=T) |
            grepl("geometric",method,ignore.case=T)|
            grepl("MSE",method,ignore.case=T)){
            TopFit<-order(UFit,decreasing=F)[seq_len(child)]
        }else{
            TopFit<-order(UFit,decreasing=T)[seq_len(child)]
        }

    #return(population[TopFit])
    return(TopFit)
}

singleRun <- function(indiv,DNAAffinity,
    genomicProfiles,DNASequenceSet,
    ChIPScore,fitness="all"){

    if(length(indiv) ==1 & class(indiv)=="list"){
      indiv <- indiv[[1]]
    }else {
      stop("Oops somthing went wrong. Not sure what your indiv object is")
    }
        lambdaPWM(genomicProfiles)<-unlist(indiv[
            grepl("ScalingFactor",names(indiv),ignore.case=T)|
            grepl("Scaling Factor",names(indiv),ignore.case=T)|
            grepl("lambda",names(indiv),ignore.case=T)])
        # Setting PWMThreshold from GA population
        PWMThreshold(genomicProfiles)<-unlist(indiv[
            grepl("PWMThreshold",names(indiv),ignore.case=T)|
            grepl("PWM Threshhold",names(indiv),ignore.case=T)])
        # Setting Number of bound molecules from GA population

        boundMolecules(genomicProfiles)<- unlist(indiv[
            grepl("boundMolecules",names(indiv),ignore.case=T)|
            grepl("Bound molecules",names(indiv),ignore.case=T)|
            grepl("^N$",names(indiv),ignore.case=T)])

        ### compitng prediction from GA population
        ## genome wide scoring
        genomeWide <- computeGenomeWideScores(genomicProfiles,DNASequenceSet,
            DNAAffinity,verbose=FALSE)

        ## compute PWMscores
        pwmScores <- computePWMScore(genomeWide,DNASequenceSet,
            ChIPScore,DNAAffinity, verbose=FALSE)

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


      return(list("occupancy"=occup,"ChIP"=chip,"gof"=GoF))
}


setChromatinStates <- function(population,chromatinStates){
    ## probably could do with some validity checks but thats for later

    ## Let's assume that chromatin states come in the form as GRanges
    ## same form as you ahve atm

    # Right for now this will just set the built affinities to a GRanges
    affinities<- lapply(population,function(x){x[grep(pattern="cs",names(x),ignore.case=T)]})

    # selecting score in chromatin states GR
    if(length(grep(x=chromatinStates$name,pattern="cs",ignore.case=T))<1){
        stateNames<-
        chromatinStates$name <- paste0("CS",chromatinStates$name)
    }
    #### NOTE this is a problem! your CS that you are adding wont always match
    #### NEED to find a better way of assinging names. A more robust a universal way of doing it
    #### Run code as is for now BUT IT NEEDS TO BE CHANGED
    chromatinStates<-chromatinStates[which(chromatinStates$name %in% names(affinities[[1]]))]
    # Removing levels
    chromatinStates<-GRanges(seqnames=as.character(seqnames(chromatinStates)),
        ranges=IRanges(start(chromatinStates), end(chromatinStates)),
        stateID=chromatinStates$name,
        DNAaffinity=rep(0,length(chromatinStates)))

    # reorder
    states <- match(names(affinities[[1]]),unique(chromatinStates$stateID))
    states <- names(affinities[[1]])[!is.na(states)]

    # replae affinity score in CS granges object
    GRStates<-vector("list", length(affinities))
    for(i in seq_along(GRStates)){
        GRStates[[i]]<-chromatinStates
        for(j in states){
            GRStates[[i]]$DNAaffinity[which(GRStates[[i]]$stateID == j)] <- affinities[[i]][[which(names(affinities[[i]])==j)]]
        }

    }

    return(GRStates)


}


## Generating a template for starting population
## with population being the number of individuals
## parameters the parameters you are interested in using
generateStartingPopulation<-function(population,parameters,names=NULL){
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

    # Then lets sample one of these one of these pre defined parameter
    # Note that there will be somme issue if you use random values
    pop<-.generateRandomGASolution(pop)

    #return final population
    return(pop)
}
