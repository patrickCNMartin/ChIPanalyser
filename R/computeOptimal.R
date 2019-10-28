computeOptimal <- function(genomicProfiles,DNASequenceSet,ChIPScore, chromatinState = NULL,
    parameterOptions = NULL,optimalMethod = "all",rank=FALSE,returnAll=TRUE,
    peakMethod="moving_kernel",cores=1){
    #validity checking


    if(!.is.genomicProfiles(genomicProfiles)){
    stop(paste0(deparse(substitute(genomicProfiles)),
    " is not a Genomic Profile Paramter object."))
    }

    if(!is.null(parameterOptions)){
        genomicProfiles<-.updateGenomicProfiles(genomicProfiles,parameterOptions)
    }

    ## Need to set a check for optimal method otherwise shit will fuck up
    if(!any(optimalMethod %in% c("pearson","spearman","kendall","ks","geometric","fscore","all"))){
        stop("Optimal method should be one of the following :
        pearson,spearman,kendall,ks,geometric,fscore  or all")
    }

    #Setting Parameters for computation
    if(all(BPFrequency(genomicProfiles)==0.25)){
    BPFrequency(genomicProfiles)<-.computeBPFrequency(
        DNASequenceSet)
    }


    if(!is.null(ChIPScore)){

           loci<-loci(ChIPScore)

    } else {
      ## change that message but essentially for compute opimntal you need to used that thing
        stop(paste("Please provide a set of loci to be analysed \n",
                   "If you have used processingChIP, you can parse that object \n"))
    }

    peakMethod <- peakMethod
    #Setting default Paramters if not provided by user
    if(length(lambdaPWM(genomicProfiles))<2){
    lambdaPWM(genomicProfiles) <- seq(0.25,5,by=0.25)
    }
    if(length(boundMolecules(genomicProfiles))<2){
    boundMolecules(genomicProfiles) <- c(1, 10, 20, 50, 100,
        200, 500,1000,2000, 5000,10000,20000,50000, 100000,
        200000, 500000, 1000000)
    }

    ## Changing cores if necessary

    if(length(loci)<cores){
        cores<-length(loci)
        warning("Number of cores requested higher than number of loci provided - some cores will be dropped")
      }

    message("Computing Genome Wide PWM Score")
    gw <- computeGenomeWideScores(genomicProfiles=genomicProfiles,
        DNASequenceSet = DNASequenceSet,
        chromatinState = chromatinState,cores=cores,
        verbose = FALSE)

    message("Computing PWM Score at Loci & Extracting Sites Above Threshold")
    pwm <- computePWMScore(DNASequenceSet = DNASequenceSet,
        genomicProfiles = gw, loci = loci,
        chromatinState = chromatinState,cores=cores, verbose = FALSE)

    message("Computing Occupancy")
    occup <- computeOccupancy(genomicProfiles = pwm,verbose = FALSE)
    # purge unused variable
   .cleanUpAfterYourself(gw,pwm)

    message("Computing ChIP-seq-like Profile")
    chip <- computeChIPProfile(loci = loci,
        genomicProfiles = occup,
        parameterOptions = parameterOptions,
        norm = TRUE ,method=peakMethod, peakSignificantThreshold= NULL,
        cores=cores,verbose = FALSE)

    message("Computing Accuracy of Profile")
    accu <- profileAccuracyEstimate(ChIPScore = ChIPScore,
        genomicProfiles = chip,
        parameterOptions = parameterOptions,method="all",cores=cores)

    #browser()
    #Extracting Optimal matrix from Profile AccuracyEstimate

    opti<-.optimalExtraction(accu,rank=FALSE)
    if(returnAll){
       return(list("Optimal"=opti,"Occupancy"=occup,"ChIPProfiles"=chip,"goodnessOfFit"=accu))
    }else{
       return(opti)
    }



}
