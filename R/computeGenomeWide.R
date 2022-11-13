################################################################################
##################### Compute genome Wide scores and metrics ###################
################################################################################


computeGenomeWideScores <- function(genomicProfiles,DNASequenceSet,
                                   chromatinState=NULL,parameterOptions=NULL,
                                   cores=1, verbose=TRUE)
  {
    # Validity Checking for each argument
    if(!is(DNASequenceSet,"BSgenome") &
    !is(DNASequenceSet,"DNAStringSet")){
    stop(paste0(deparse(substitute(DNASequenceSet)),
    " is not a BSgenome Object or a DNAStringSet"))
    }

    if(is(DNASequenceSet,"BSgenome")){
    DNASequenceSet<-getSeq(DNASequenceSet)
    }

    if(!.is.genomicProfiles(genomicProfiles)){
    stop(paste0(deparse(substitute(genomicProfiles)),
    " is not a Genomic Profile Paramter object.
    Please set Genomic Profile Parameters."))
    }

    if(!is(chromatinState, "GRanges")){
    stop(paste0(chromatinState," must be a GRanges Object."))
    }

    if(!is.null(parameterOptions)){
        genomicProfiles<-.updateGenomicProfiles(genomicProfiles,parameterOptions)
    }
    # Extracting genomic Profile Parameters
    lambda <- lambdaPWM(genomicProfiles)
    PWMMat <- PositionWeightMatrix(genomicProfiles)
    strand <- whichstrand(genomicProfiles)
    strandRule <- strandRule(genomicProfiles)
    # Extracting DNA Accessibility in parallel (still dont know why it chrashes with hg19)
    if(!is.null(chromatinState)){
        ## Factor correction  and ID match

        DNASequenceSet<-DNASequenceSet[which(names(DNASequenceSet) %in% as.character(seqnames(chromatinState)))]
        ## Split by Chromosome
        DNA<-split(DNASequenceSet,names(DNASequenceSet))
        Access<-split(chromatinState, seqnames(chromatinState))

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
    message("Extracting genome wide scores \n")
        if(!is.null(chromatinState)){
        if(strand == "+") message("Considering Chromatin State ~ positive strand \n")
        if(strand == "-") message("Considering Chromatin State ~ negative strand \n")
        if(strand == "+-" | strand == "-+"){
        message("Considering Chromatin State ~ Both strands \n")
        }
      } else {
        if(strand == "+") message("Whole Genome ~ positive strand \n")
        if(strand == "-") message("Whole Genome ~ negative strand \n")
        if(strand == "+-" | strand == "-+"){
        message("Whole Genome ~ Both strands \n")
        }
      }
    }


    ### Custom DNSStringSet Split based on number of cores

    buffer <- .splitRanges(DNASequenceSet,cores)
    DNASequenceSet<-buffer$rangeSet
    cores<-buffer$cores

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
    .maxPWMScore(genomicProfiles)<-maxPWMScore
    .minPWMScore(genomicProfiles)<-minPWMScore
    .averageExpPWMScore(genomicProfiles)<-averageExpPWMScore
    .DNASequenceLength(genomicProfiles)<-DNASequenceLength
    if(!is.null(chromatinState)){
    .tags(genomicProfiles)<-"genomeWideCS"
    .paramTag(genomicProfiles)<-"genomeWide"
  } else{
    .tags(genomicProfiles)<-"genomeWide"
    .paramTag(genomicProfiles)<-"genomeWide"
  }
    #genomicProfiles <-.maxPWMScoreReplace(genomicProfileParameters,
      #  maxPWMScore)
    #genomicProfiles <-.minPWMScoreReplace(genomicProfileParameters,
        #minPWMScore)
    #genomicProfiles <-.averageExpPWMScoreReplace(
        #genomicProfileParameters,averageExpPWMScore)
  #  DNASequenceLength(genomicProfileParameters) <- DNASequenceLength

    return(genomicProfiles)

}
