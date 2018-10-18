#######################################################################
# Functions
#######################################################################
computeGenomeWidePWMScore<-function(DNASequenceSet,
    genomicProfileParameters, DNAAccessibility = NULL,cores=1,
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
    #if(!is.null(DNAAccessibility)){
    #DNASequenceSet <- getSeq(DNASequenceSet,DNAAccessibility)
    #names(DNASequenceSet) <- seqnames(DNAAccessibility)
    #DNASequenceSet <- DNASequenceSet[which(sapply(DNASequenceSet,length) >
    #    ncol(PWMMat))]
    #}

    # Extracting DNA Accessibility in parallel (still dont know why it chrashes with hg19)
    if(!is.null(DNAAccessibility)){
        ## Factor correction  and ID match

        DNASequenceSet<-DNASequenceSet[which(names(DNASequenceSet) %in% as.character(seqnames(DNAAccessibility)))]
        ## Split by Chromosome
        DNA<-split(DNASequenceSet,names(DNASequenceSet))
        Access<-split(DNAAccessibility, seqnames(DNAAccessibility))

        ## Matching seqlevels
        Access<- Access[match(names(DNA), names(Access))]


        ## parallel Intersect for large genomes
        DNASequenceSet<-parallel::mcmapply(.internalDNASequenceFromAccess,DNA=DNA,Access=Access,mc.cores=cores)


        # Rebuilding DNASequenceSet
        if(length(DNASequenceSet)>1){
            buffer<-DNASequenceSet[[1]]
            for(sub in 2:length(DNASequenceSet)){
                buffer<-c(buffer,DNASequenceSet[[sub]])
            }
            DNASequenceSet<-buffer
        } else {
            DNASequenceSet<-DNASequenceSet[[1]]
        }

        DNASequenceSet <- DNASequenceSet[which(width(DNASequenceSet) > ncol(PWMMat))]

    }
    # Computing DNA sequance Length
    DNASequenceLength <- sum(as.numeric(width(DNASequenceSet)))

    # Progress Messages when required.
    if(verbose){
    message("Scoring whole genome \n")
        if(strand == "+") message("Accessible DNA ~ positive strand \n")
        if(strand == "-") message("Accessible DNA ~ negative strand \n")
        if(strand == "+-" | strand == "-+"){
        message("Accessible DNA ~ Both strands \n")
        }
    }


    ### Custom DNSStringSet Split based on number of cores

    DNASequenceSet <- .splitDNARanges(DNASequenceSet,cores)

    ### Computing PWM Score parallel
    DNASequenceScoreSetTotalAcces <- parallel::mclapply(DNASequenceSet,
      .scoreDNAStringSet,PWM=PWMMat,
       strand=strand,strandRule=strandRule,mc.cores=cores)


    DNASequenceScoreSetTotalAcces <- unlist(lapply(DNASequenceScoreSetTotalAcces,"[[",1))

    #Message printing when required
    if(verbose){
    message("Computing Mean waiting time \n")
    }

    # Compute mean waiting time, max PWM score and min PWM score

    averageExpPWMScore <- rep(0,length(lambda))
    sumExpPWMScoreLocal <- rep(0,length(lambda))

    ## This is also a limiting part. parallel as well?
    ## Clean your parallel R script and comment it

    for(i in seq_along(lambda)){
        sumExpPWMScoreLocal[i] <- sum(exp(DNASequenceScoreSetTotalAcces * (1/lambda[i])))
        averageExpPWMScore[i] <- sumExpPWMScoreLocal[i]/DNASequenceLength
    }

    maxPWMScore <- max(DNASequenceScoreSetTotalAcces)
    minPWMScore <- min(DNASequenceScoreSetTotalAcces)

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
