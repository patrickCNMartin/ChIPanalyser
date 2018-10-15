#######################################################################
# Functions
#######################################################################

computePWMScore <- function(DNASequenceSet,genomicProfileParameters,
    setSequence = NULL,DNAAccessibility = NULL,cores=1,verbose = TRUE){

    # Validity checking
    if(!.is.genomicProfileParameter(genomicProfileParameters)){
    stop(paste0(deparse(substitute(genomicProfileParameters)),
    " is not a Genomic Profile Paramter object.
    Please set a Genomic Profile Parameters Object."))
    }
    if(!.is.genomeWideComputed(genomicProfileParameters)){
    stop(paste0("Genome Wide PWM Score has not yet been computed in ",
    deparse(substitute(genomicProfileParameters))))
    }
    if(class(DNASequenceSet) != "BSgenome" &
    class(DNASequenceSet)!="DNAStringSet"){
    stop(paste0(deparse(substitute(DNASequenceSet)),
    " is not a BSgenome Object or a DNAStringSet"))
    }
    if(class(DNASequenceSet) == "BSgenome"){
    DNASequenceSet<-getSeq(DNASequenceSet)
    }
    if(class(DNAAccessibility) != "GRanges" & !is.null(DNAAccessibility)){
    stop(paste0(deparse(substitute(DNAAccessibility)),
    " must be a GRanges Object."))
    }
    if(class(setSequence) != "GRanges" & !is.null(setSequence)){
    stop(paste0(deparse(substitute(setSequence)),
    " must be a GRanges Object."))
    }

    #Extracting genomic Profile Parameters
    PWMThreshold <- PWMThreshold(genomicProfileParameters)
    minPWMScore <- minPWMScore(genomicProfileParameters)
    maxPWMScore <- maxPWMScore(genomicProfileParameters)
    PWM <- PositionWeightMatrix(genomicProfileParameters)
    strand <- whichstrand(genomicProfileParameters)
    strandRule <- strandRule(genomicProfileParameters)


    #Calculating Threshhold Value
    PWMThresholdLocal <- minPWMScore + PWMThreshold*(maxPWMScore-minPWMScore)

    #Processing DNAAccessibility
    if(!is.null(DNAAccessibility)){
        if(length(DNAAccessibility$DNAAccessibility)<1){
            DNAAccessibility$DNAAccessibility <- rep(1, length(DNAAccessibility))
        }
    }

    ### setSequence split

    if(length(setSequence)>cores & cores>1){
        message("Multi Core Sequence Split")
        setSequence<-setSequence[which(width(setSequence)>2*ncol(PWM))]
        localSequence <- .splitDNARanges(setSequence,cores)
        message("Multi Core PWM Scores Extraction")
        ## Computing PWM Score above threshold
        Scores<-parallel::mclapply(localSequence,.internalPWMScoreExt,
            DNASequenceSet=DNASequenceSet,
            DNAAccessibility=DNAAccessibility,
            PWM=PWM,PWMThresholdLocal=PWMThresholdLocal,minPWMScore=minPWMScore,
            maxPWMScore=maxPWMScore,strand=strand,strandRule=strandRule,mc.cores=cores)

        AccessibleSequence<-unlist(lapply(Scores,"[[",1))
        NoAcc<-unlist(lapply(Scores,"[[",2))

    } else {
        message("Single Core PWM Scores Extraction")
        localSequence<-setSequence[which(width(setSequence)>2*ncol(PWM)+1)]
        ## Computing PWM Score above threshold
        Scores<-.internalPWMScoreExt(localSequence,
            DNASequenceSet=DNASequenceSet,
            DNAAccessibility=DNAAccessibility,
            PWM=PWM,PWMThresholdLocal=PWMThresholdLocal,minPWMScore=minPWMScore,
            maxPWMScore=maxPWMScore,strand=strand,strandRule=strandRule)
        AccessibleSequence<-Scores[[1]]
        NoAcc<-Scores[[2]]
    }

    ## Re-assigning names

    names(AccessibleSequence)<-sapply(AccessibleSequence,function(x)names(x)[1])

    genomicProfileParameters <-.AllSitesAboveThresholdReplace(
        genomicProfileParameters,GRangesList(
        AccessibleSequence[which(lapply(AccessibleSequence,length) != 0)]))

    genomicProfileParameters <-.NoAccessReplace(genomicProfileParameters,
        NoAcc)
    return(genomicProfileParameters)

}
