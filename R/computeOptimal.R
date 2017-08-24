computeOptimal <- function(DNASequenceSet,genomicProfileParameters,
    LocusProfile,setSequence, DNAAccessibility = NULL,
    occupancyProfileParameters = NULL,parameter = "all",
    peakMethod="moving_kernel"){
    #validity checking


    if(!.is.genomicProfileParameter(genomicProfileParameters)){
    stop(paste0(deparse(substitute(genomicProfileParameters)),
    " is not a Genomic Profile Paramter object."))
    }

    if(!is.null(occupancyProfileParameters)){
        if (!.is.occupancyProfileParameters(occupancyProfileParameters)){
        stop(paste0(deparse(substitute(occupancyProfileParameters)),
        " is not a Genomic Profile Paramter object."))
        }
    } else{
        occupancyProfileParameters <- occupancyProfileParameters()
    }

    #Setting Parameters for computation
    if(all(BPFrequency(genomicProfileParameters)==0.25)){
    BPFrequency(genomicProfileParameters)<-.computeBPFrequency(
        DNASequenceSet)
    } else {
    BPFrequency(genomicProfileParameters)<-BPFrequency(
        genomicProfileParameters)
    }
    peakMethod <- peakMethod
    #Setting default Paramters if not provided by user
    if(length(ScalingFactorPWM(genomicProfileParameters))<2){
    ScalingFactorPWM(genomicProfileParameters) <- c(0.25, 0.5, 0.75, 1, 1.25,
        1.5, 1.75, 2, 2.5, 3, 3.5 ,4 ,4.5, 5)
    }
    if(length(boundMolecules(occupancyProfileParameters))<2){
    boundMolecules(occupancyProfileParameters) <- c(1, 10, 20, 50, 100,
        200, 500,1000,2000, 5000,10000,20000,50000, 100000,
        200000, 500000, 1000000)
    }

    message("Computing Genome Wide PWM Score")
    GenomeWide <- computeGenomeWidePWMScore(DNASequenceSet = DNASequenceSet,
        genomicProfileParameters = genomicProfileParameters,
        DNAAccessibility = DNAAccessibility,
        verbose = FALSE)

    message("Computing PWM Score at Loci & Extracting Sites Above Threshold")
    SiteSpecific <- computePWMScore(DNASequenceSet = DNASequenceSet,
        genomicProfileParameters = GenomeWide, setSequence = setSequence,
        DNAAccessibility = DNAAccessibility, verbose = FALSE)

    message("Computing Occupancy")
    Occupancy <- computeOccupancy(AllSitesPWMScore = SiteSpecific,
        occupancyProfileParameters = occupancyProfileParameters,
        verbose = FALSE)

    message("Computing ChIP-seq-like Profile")
    PedictedProfile <- computeChipProfile(setSequence = setSequence,
        occupancy = Occupancy,
        occupancyProfileParameters = occupancyProfileParameters,
        norm = TRUE ,method=peakMethod, peakSignificantThreshold= NULL,
        verbose = FALSE)

    message("Computing Accuracy of Profile")
    ProfileAccuracy <- profileAccuracyEstimate(LocusProfile = LocusProfile,
        predictedProfile = PedictedProfile,
        occupancyProfileParameters = occupancyProfileParameters)


    #Extracting Optimal matrix from Profile AccuracyEstimate
    AllMatrix <- vector("list",3)
    AllParam <- vector("list",3)
    Param <- c("meanCorr","meanMSE","meanTheta")
    for(i in 1:3){
        AllParam[[i]] <- c(0,0)
        AllMatrix[[i]] <- matrix(0,nrow=length(ScalingFactorPWM(
            genomicProfileParameters)),
            ncol=length(boundMolecules(occupancyProfileParameters)))
        rownames(AllMatrix[[i]]) <- ScalingFactorPWM(genomicProfileParameters)
        colnames(AllMatrix[[i]]) <- boundMolecules(occupancyProfileParameters)

        ParamSplit <- lapply(ProfileAccuracy,function(x){x<-x[[1]][Param[i]]})
        ParamSplit <- split(ParamSplit,1:length(boundMolecules(
            occupancyProfileParameters)))

        for(j in 1:length(ParamSplit)){
            AllMatrix[[i]][,j] <- unname(unlist(unname(ParamSplit[[j]])))
        }
        if(Param[i]=="meanMSE"){
            buffer <- which(AllMatrix[[i]]==min(AllMatrix[[i]]),arr.ind=TRUE)
        }else{
            buffer <- which(AllMatrix[[i]]==max(AllMatrix[[i]]),arr.ind=TRUE)
        }

        AllParam[[i]][1] <- rownames(AllMatrix[[i]])[buffer[,1]]
        AllParam[[i]][2] <- colnames(AllMatrix[[i]])[buffer[,2]]
    }

    names(AllMatrix) <- Param
    names(AllParam) <- Param
    message("Extracting Optimal Set of Parameters")
    if(parameter == "correlation"){
        optimalParam <- AllParam[["meanCorr"]]
        optimalMatrix <- AllMatrix[["meanCorr"]]
    }
    if(parameter == "MSE"){
        optimalParam <- AllParam[["meanMSE"]]
        optimalMatrix <- AllMatrix[["meanMSE"]]
    }
    if(parameter == "theta"){
        optimalParam <- AllParam[["meanTheta"]]
        optimalMatrix <- AllMatrix[["meanTheta"]]
    }
    if(parameter == "all"){
        optimalParam <- AllParam
        optimalMatrix <- AllMatrix
    }
    return(list("Optimal Parameters" = optimalParam,
    "Optimal Matrix" = optimalMatrix,"Parameter" = parameter))

}
