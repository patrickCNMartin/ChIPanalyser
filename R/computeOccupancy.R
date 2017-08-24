#######################################################################
# Functions
#######################################################################


computeOccupancy <- function(AllSitesPWMScore,
    occupancyProfileParameters=NULL,
    norm=TRUE,verbose=TRUE) {
    # Validity checking
    if(!.is.occupancyProfileParameters(occupancyProfileParameters) &
    !is.null(occupancyProfileParameters)){
    stop(paste0(deparse(substitute(occupancyProfileParameters)),
    "is not an occupancyProfileParamaters Object."))
    }
    if (!.is.genomicProfileParameter(AllSitesPWMScore)){
    stop(paste0(deparse(substitute(AllSitesPWMScore)),
    " is not a Genomic Profile Paramter object."))
    }
    if(!.is.PWMScoreComputed(AllSitesPWMScore)){
    stop(paste0("PWM Score at sites higher than threshold " ,
    "have not yet been computed in",
    deparse(substitute(AllSitesPWMScore))))
    }
    #Generating occupancyProfileParamters if not given by user.
    #All Value will be defaut settings
    if(is.null(occupancyProfileParameters)){
    occupancyProfileParameters <- occupancyProfileParameters()
    }

    # GenomicProfileParameter Parameter Extraction
    PWMGRList <- AllSitesAboveThreshold(AllSitesPWMScore)
    DNALength <- as.numeric(DNASequenceLength(AllSitesPWMScore))
    averageExpPWMScore <- averageExpPWMScore(AllSitesPWMScore)
    lambda <- ScalingFactorPWM(AllSitesPWMScore)

    #Occupancy Profile Parameter Extraction
    ploidy <- as.numeric(ploidy(occupancyProfileParameters))
    boundMolecules <- boundMolecules(occupancyProfileParameters)
    backgroundSignal <- backgroundSignal(occupancyProfileParameters)
    maxSignal <- maxSignal(occupancyProfileParameters)



    # Computing Occupancy at sites higher than threshold
    MultiParam <- vector("list", (length(lambda)*length(boundMolecules)))
    PWMScore <- vector("list", length(PWMGRList))
    # Extracting names of regions
    name <- c()
    for(i in 1:length(PWMGRList)){
        name <- c(name,rep(names(PWMGRList)[[i]],length(PWMGRList[[i]])))
    }

    Occupancy <- vector("list",length(PWMGRList))
    names(Occupancy) <- names(PWMGRList)
    result <- GRangesList()
    emptyGR <- GRanges()
    # Progress message when required
    if(verbose){
    message("Computing Occupancy at sites higher than threshold. \n")
    }

    counter <- 0
    ParaVal <- c()

    #Computing Occupancy
    for(k in 1:length(lambda)){
        for(j in 1:length(boundMolecules)){
            PWMScore <- as.vector(as.matrix(mcols(unlist(PWMGRList))))
            Occupancy <- rep(0,length(PWMScore))
            Occupancy <- (boundMolecules[j]*exp((1/lambda[k])*PWMScore))/
                (boundMolecules[j]*exp((1/lambda[k])*PWMScore) +
                DNALength*ploidy*averageExpPWMScore[k])

            Occupancy <- backgroundSignal + Occupancy*
                (maxSignal-backgroundSignal)
    # Normalising Ocupancy Signal
            if(norm == TRUE){
                maxOccupancy <- max(c(maxSignal,max(Occupancy)))
                Occupancy <- Occupancy/maxOccupancy
                ZeroBackground <- backgroundSignal/maxOccupancy
            } else {
                ZeroBackground <- backgroundSignal
            }

            buffer <- unlist(PWMGRList)
            buffer$Occupancy <- Occupancy
    # Extracting and pasting names of Parameters
            names(buffer) <- name
            if(length(name > 1)){
                level <- factor(names(buffer),levels=unique(name))
                result <- split(buffer,level)
            } else {
                result <-  buffer
            }

            counter <- counter+1
            MultiParam[[counter]] <- result
            ParaVal <- c(ParaVal,paste0("lambda = ",lambda[k]," & ",
            "boundMolecules = ",boundMolecules[j]))
        }
    }
    # ReBuilding GenomicProfileParameters
    names(MultiParam) <- ParaVal

    AllSitesPWMScore<-.ZeroBackgroundReplace(AllSitesPWMScore,
        ZeroBackground)
    AllSitesPWMScore<-.AllSitesAboveThresholdReplace(AllSitesPWMScore,
        MultiParam)

    return(AllSitesPWMScore)
}
