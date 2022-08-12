#######################################################################
# Functions
#######################################################################

computePWMScore <- function(genomicProfiles,DNASequenceSet,
    loci = NULL,chromatinState = NULL,parameterOptions=NULL,cores=1,verbose = TRUE){

    # Validity checking
    # this needs to be re done but low priority at this point

    if(!.is.genomicProfiles(genomicProfiles)){
    stop(paste0(deparse(substitute(genomicProfiles)),
    " is not a Genomic Profile Paramter object.
    Please set a Genomic Profile Parameters Object."))
    }

    if(class(DNASequenceSet) != "BSgenome" &
    class(DNASequenceSet)!="DNAStringSet"){
    stop(paste0(deparse(substitute(DNASequenceSet)),
    " is not a BSgenome Object or a DNAStringSet"))
    }
    if(class(DNASequenceSet) == "BSgenome"){
    DNASequenceSet<-getSeq(DNASequenceSet)
    }
    if(class(chromatinState) != "GRanges" & !is.null(chromatinState)){
    stop(paste0(deparse(substitute(chromatinState)),
    " must be a GRanges Object."))
    }

    if(!is.null(parameterOptions)){
        genomicProfiles<-.updateGenomicProfiles(genomicProfiles,parameterOptions)
    }


    # Loci check from chipscore
    if(!is.null(loci)){
        if(class(loci)=="ChIPScore"){
           loci<-loci(loci)
        } else if(class(loci)=="Granges") {
           loci <-loci
        }
    } else {
        stop(paste("Please provide a set of loci to be analysed \n",
                   "If you have used processingChIP, you can parse that object \n",
                   "Otherwise, you can provide your own loci as a GRanges. "))
    }
    if(length(loci)<cores){
      cores<-length(loci)
      warning("Number of cores requested higher than number of loci provided - some cores will be dropped")
    }
    #Extracting genomic Profile Parameters
    PWMThreshold <- PWMThreshold(genomicProfiles)
    minPWMScore <- minPWMScore(genomicProfiles)
    maxPWMScore <- maxPWMScore(genomicProfiles)
    PWM <- PositionWeightMatrix(genomicProfiles)
    strand <- whichstrand(genomicProfiles)
    strandRule <- strandRule(genomicProfiles)


    #Calculating Threshhold Value
    PWMThresholdLocal <- minPWMScore + PWMThreshold*(maxPWMScore-minPWMScore)

    #Processing chromatinState
    if(!is.null(chromatinState)){
        if(length(chromatinState$DNAaffinity)<1){
            chromatinState$DNAaffinity <- rep(1, length(chromatinState))
        }
    }

    ### loci split

    if(length(loci)>=cores & cores>1){
        #message("Multi Core Sequence Split")
        loci<-loci[which(width(loci)>2*ncol(PWM))]
        buffer <- .splitRanges(loci,cores)
        localSequence<-buffer$rangeSet
        cores<-buffer$cores
        if(verbose){
            message("PWM Scores Extraction")
        }
        
        ## Computing PWM Score above threshold
        Scores<-parallel::mclapply(localSequence,.internalPWMScoreExt,
            DNASequenceSet=DNASequenceSet,
            chromatinState=chromatinState,
            PWM=PWM,PWMThresholdLocal=PWMThresholdLocal,minPWMScore=minPWMScore,
            maxPWMScore=maxPWMScore,strand=strand,strandRule=strandRule,mc.cores=cores)

        AccessibleSequence<-Scores


    } else {
        if(verbose){
            message("PWM Scores Extraction")
        }
        localSequence<-loci[which(width(loci)>2*ncol(PWM)+1)]
        ## Computing PWM Score above threshold
        Scores<-.internalPWMScoreExt(localSequence,
            DNASequenceSet=DNASequenceSet,
            chromatinState=chromatinState,
            PWM=PWM,PWMThresholdLocal=PWMThresholdLocal,minPWMScore=minPWMScore,
            maxPWMScore=maxPWMScore,strand=strand,strandRule=strandRule)
        AccessibleSequence<-Scores


    }

    ## Re-assigning names
     AccessibleSequence<-unlist(AccessibleSequence)
    .profiles(genomicProfiles)<-AccessibleSequence

     lociDrop<-names(AccessibleSequence)[which(sapply(AccessibleSequence,length)==0)]
     if(length(lociDrop)!=0){
        .drop(genomicProfiles) <- lociDrop
     } else {
        .drop(genomicProfiles) <-"No loci dropped"
     }
     .tags(genomicProfiles)<-"PWMScore"
     .paramTag(genomicProfiles)<-"PWMScore"

    return(genomicProfiles)

}
