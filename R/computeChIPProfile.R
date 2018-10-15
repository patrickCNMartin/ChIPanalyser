#######################################################################
# Functions
#######################################################################


computeChipProfile <- function(setSequence,occupancy,
    occupancyProfileParameters = NULL,
    norm = TRUE, method = c("moving_kernel","truncated_kernel","exact"),
    peakSignificantThreshold = NULL,cores=1, verbose = TRUE){

    # Validity checking
    if(class(setSequence)!="GRanges"){
    stop(paste0(deparse(substitute(SetSequence)),
    " must be a GRanges Object"))
    }
    if(!.is.genomicProfileParameter(occupancy)){
    stop(paste0(deparse(substitute(occupancy)),
    " is not a Genomic Profile Paramter Object"))
    }
    if(!.is.occupancyProfileParameters(occupancyProfileParameters) &
    !is.null(occupancyProfileParameters)){
    stop(paste0(deparse(substitute(OPP)),
    " is not an occupancyProfileParamaters Object."))
    }
    if(is.null(occupancyProfileParameters)){
    occupancyProfileParameters <- occupancyProfileParameters()
    }
    if(cores > detectCores()){
        message(paste0(detectCores(), " cores are available on this machine.
        Tying to set higher number of cores than available"))
    }
    # Extraction of Ocuupancy and associated values
    Occup <- AllSitesAboveThreshold(occupancy)
    ZeroBackground <- .ZeroBackground(occupancy)

    #Extraction of Occupancy Profile Parameters
    stepSize <- stepSize(occupancyProfileParameters)
    backgroundSignal <- backgroundSignal(occupancyProfileParameters)
    removeBackground <- removeBackground(occupancyProfileParameters)
    chipMean <- chipMean(occupancyProfileParameters)
    chipSd <- chipSd(occupancyProfileParameters)

    if(!is.null(chipSmooth(occupancyProfileParameters))){
        chipSmooth <- chipSmooth(occupancyProfileParameters)
    } else {
        chipSmooth <- NULL
    }
    maxSignal <- maxSignal(occupancyProfileParameters)


    # Extracting names of sequences with no accesible DNA
    NoAccess <- names(setSequence)[(names(setSequence) %in%
        names(Occup[[1]])==FALSE)]
    if(length(NoAccess) > 0){
        cat("No Profile for:",NoAccess,
        "  --  Do Not Contain Accessible Sites","\n", sep=" ")
    }

    # SetSequence fragmentation
    #Spliting setSequence for Chip PRofile computing
    setSequence <- setSequence[names(setSequence) %in% names(Occup[[1]])]

    LocalSet <- split(setSequence, seq_along(setSequence))

    #Computing Chip like profile
    if(verbose){
        message("Computing ChIP Profile \n")
    }

    SplitGRList <-parallel::mclapply(LocalSet,.internalChIPLociSplit,stepSize,mc.cores=cores)

    names(SplitGRList) <- names(setSequence)

    OccupancyVals <- lapply(Occup,function(x){x<-lapply(x,function(y){
        y<-as.numeric(as.matrix(mcols(y)[,3]))})})

    ##method set
    if(length(Occup)>length(Occup[[1]])){
        OccupancyVals <- parallel::mclapply(Occup,.internalChIPOccupValsParam,mc.cores=cores)
        profile <- parallel::mcmapply(.internalChIPParam,Occup=Occup,
        OccupancyVals=OccupancyVals,
        MoreArgs=list(SplitGRList=SplitGRList,LocalSet=LocalSet,
        chipMean=chipMean,chipSd=chipSd,
        stepSize=stepSize,norm=norm,chipSmooth=chipSmooth,
        peakSignificantThreshold=peakSignificantThreshold,
        ZeroBackground=ZeroBackground,
        removeBackground=removeBackground,method=method),
        mc.cores=cores,SIMPLIFY=FALSE)
    } else {
    profile <- vector("list", length(Occup))
    OccupancyVals <- vector("list", length(Occup))
    for(i in seq_along(Occup)){
        OccupancyVals[[i]]<- parallel::mclapply(Occup[[i]],.internalChIPOccupValsLoci,mc.cores=cores)
        profile[[i]]<-SplitGRList
        profile[[i]] <- parallel::mcmapply(.internalChIPLoci,
        profile=profile[[i]],Occup=Occup[[i]],LocalSet=LocalSet,
        OccupancyVals=OccupancyVals[[i]],
        MoreArgs=list(chipMean=chipMean,chipSd=chipSd,
        stepSize=stepSize,norm=norm,chipSmooth=chipSmooth,
        peakSignificantThreshold=peakSignificantThreshold,
        ZeroBackground=ZeroBackground,
        removeBackground=removeBackground,method=method), mc.cores=cores)

    }
    }

    names(profile)<-names(Occup)
    return(profile)

}
