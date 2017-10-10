#######################################################################
# Functions
#######################################################################


computeChipProfile <- function(setSequence,occupancy,
    occupancyProfileParameters = NULL,
    norm = TRUE, method = c("moving_kernel","truncated_kernel","exact"),
    peakSignificantThreshold = NULL, verbose = TRUE){

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
    SplitGRList <- vector("list",length(setSequence))
    names(SplitGRList) <- names(setSequence)
    for(i in seq_along(SplitGRList)){
        stepIndex <- seq(from=1, to=width(setSequence[i]), by=stepSize)
        SplitSeq <- GRanges(seqnames=S4Vectors::Rle(as.character(
            seqnames(setSequence[i])),
            length(stepIndex)),
            ranges = IRanges::IRanges(start = start(setSequence[i])+
                stepIndex-1,end =start(setSequence[i])+stepIndex-1+(stepSize)),
            strand = "*", ChIP=rep(0,length(stepIndex)))

        names(SplitSeq) <- rep(names(setSequence[i]), length(SplitSeq))
        SplitGRList[[i]]<- SplitSeq
    }
    #Computing Chip like profile
    if(verbose){
        message("Computing ChIP Profile \n")
    }
    profile <- vector("list",length(Occup))

    OccupancyVals <- lapply(Occup,function(x){x<-lapply(x,function(y){
        y<-as.numeric(as.matrix(mcols(y)[,2]))})})

    for(i in seq_along(Occup)){
        profile[[i]]<-SplitGRList

        for(j in seq_along(Occup[[i]])){
            stepIndex<-seq(from=1, to=width(setSequence[j]), by=stepSize)
            occupancyAbundanceChIPLocal<-rep(
                ZeroBackground,width(setSequence[j]))
            occupancyAbundanceChIPLocal[(start(Occup[[i]][[j]]) -
                start(setSequence[j]) + 1)] <- OccupancyVals[[i]][[j]]
            occupancyAbundanceChIPLocal[which(occupancyAbundanceChIPLocal <
                removeBackground)] <- 0
            if(any(method=="exact")){
            occupancyAbundanceChIPLocal <- .generateChIPProfile(
                occupancyAbundanceChIPLocal,chipMean,chipSd,
                chipSmooth,norm=norm, quick=FALSE,
                peakSignificantThreshold=peakSignificantThreshold)
            }else if(any(method=="truncated_kernel")){
            occupancyAbundanceChIPLocal <- .generateChIPProfile(
                occupancyAbundanceChIPLocal,chipMean,chipSd,
                chipSmooth,norm=norm, quick=TRUE,
                peakSignificantThreshold=peakSignificantThreshold)
            } else if(any(method=="moving_kernel")){
            occupancyAbundanceChIPLocal <- .generateChIPProfileRcpp(
                occupancyAbundanceChIPLocal, chipMean, chipSd, chipSmooth,
                norm = norm,
                peakSignificantThreshold=peakSignificantThreshold)
            }
            profile[[i]][[j]]$ChIP<-occupancyAbundanceChIPLocal[stepIndex]
        }
    }
    names(profile)<-names(Occup)
    return(profile)

}
