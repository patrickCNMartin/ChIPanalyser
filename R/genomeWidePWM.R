#######################################################################
# Functions
#######################################################################
computeGenomeWidePWMScore<-function(DNASequenceSet,
    genomicProfileParameters, DNAAccessibility = NULL,
    verbose = TRUE){

    # Validity Checking for each argument
    if(class(DNASequenceSet) != "BSgenome" &
    class(DNASequenceSet) != "DNAStringSet"){
    stop(paste0(deparse(substitute(DNASequenceSet)),
    " is not a BSgenome Object or a DNAStringSet"))
    }

    if(class(DNASequenceSet) == "BSgenome"){
    DNASequenceSet<-getSeq(DNASequenceSet)
    }

    if(!.is.genomicProfileParameter(genomicProfileParameters)){
    stop(paste0(deparse(substitute(genomicProfileParameters)),
    " is not a Genomic Profile Paramter object.
    Please set Genomic Profile Parameters."))
    }

    if(class(DNAAccessibility) != "GRanges" & !is.null(DNAAccessibility)){
    stop(paste0(DNAAccessibility," must be a GRanges Object."))
    }

    # Extracting genomic Profile Parameters
    lambda <- ScalingFactorPWM(genomicProfileParameters)
    PWMMat <- PositionWeightMatrix(genomicProfileParameters)
    strand <- whichstrand(genomicProfileParameters)
    strandRule <- strandRule(genomicProfileParameters)

    # Extracting DNA Accessibility
    if(!is.null(DNAAccessibility)){
    DNASequenceSet <- getSeq(DNASequenceSet,DNAAccessibility)
    names(DNASequenceSet) <- seqnames(DNAAccessibility)
    DNASequenceSet <- DNASequenceSet[which(sapply(DNASequenceSet,length) >
        ncol(PWMMat))]
    }
    # Progress Messages when required.
    if(verbose){
    message("Scoring whole genome \n")
        if(strand == "+") message("Accessible DNA ~ positive strand \n")
        if(strand == "-") message("Accessible DNA ~ negative strand \n")
        if(strand == "+-" | strand == "-+"){
        message("Accessible DNA ~ Both strands \n")
        }
    }

    # Computing PWM Score
    DNASequenceScoreSetTotalAcces <- .scoreDNAStringSet(PWMMat,
        DNASequenceSet,strand = strand,strandRule = strandRule)
    DNASequenceScoreSetTotalAcces <- DNASequenceScoreSetTotalAcces[[1]]

    #Message printing when required
    if(verbose){
    message("Computing Mean waiting time\n")
    }

    # Compute mean waiting time, max PWM score and min PWM score
    # Computing DNA sequance Length
    DNASequenceLength <- sum(unlist(
        lapply(DNASequenceScoreSetTotalAcces,length)))
    averageExpPWMScore <- rep(0,length(lambda))
    sumExpPWMScoreLocal <- vector("list", length(lambda))

    for(i in seq_along(lambda)){
        sumExpPWMScoreLocal[[i]] <- sapply(lapply(
            lapply(DNASequenceScoreSetTotalAcces,"*", (1/lambda[i])),exp),sum)
        averageExpPWMScore[i] <- sum(
            sumExpPWMScoreLocal[[i]])/DNASequenceLength
    }

    maxPWMScore <- max(unlist(lapply(DNASequenceScoreSetTotalAcces,max)))
    minPWMScore <- min(unlist(lapply(DNASequenceScoreSetTotalAcces,min)))

    #Updating GenomicProfileParameters object
    genomicProfileParameters <-.maxPWMScoreReplace(genomicProfileParameters,
        maxPWMScore)
    genomicProfileParameters <-.minPWMScoreReplace(genomicProfileParameters,
        minPWMScore)
    genomicProfileParameters <-.averageExpPWMScoreReplace(
        genomicProfileParameters,averageExpPWMScore)
    DNASequenceLength(genomicProfileParameters) <- DNASequenceLength

    return(genomicProfileParameters)
}
