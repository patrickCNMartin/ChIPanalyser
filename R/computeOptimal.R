computeOptimal <- function(DNASequenceSet,genomicProfileParameters,
    LocusProfile,setSequence, DNAAccessibility = NULL,
    occupancyProfileParameters = NULL,optimalMethod = "all",
    peakMethod="moving_kernel",cores=1){
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

    ## Need to set a check for optimal method otherwise shit will fuck up
    if(!any(optimalMethod %in% c("pearson","spearman","kendall","ks","geometric","fscore","all"))){
        stop("Optimal method should be one of the following :
        pearson,spearman,kendall,ks,geometric,fscore  or all")
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
        DNAAccessibility = DNAAccessibility,cores=cores,
        verbose = FALSE)

    message("Computing PWM Score at Loci & Extracting Sites Above Threshold")
    SiteSpecific <- computePWMScore(DNASequenceSet = DNASequenceSet,
        genomicProfileParameters = GenomeWide, setSequence = setSequence,
        DNAAccessibility = DNAAccessibility,cores=cores, verbose = FALSE)

    message("Computing Occupancy")
    Occupancy <- computeOccupancy(AllSitesPWMScore = SiteSpecific,
        occupancyProfileParameters = occupancyProfileParameters,
        verbose = FALSE)

    message("Computing ChIP-seq-like Profile")
    PredictedProfile <- computeChipProfile(setSequence = setSequence,
        occupancy = Occupancy,
        occupancyProfileParameters = occupancyProfileParameters,
        norm = TRUE ,method=peakMethod, peakSignificantThreshold= NULL,
        cores=cores,verbose = FALSE)

    message("Computing Accuracy of Profile")
    ProfileAccuracy <- profileAccuracyEstimate(LocusProfile = LocusProfile,
        predictedProfile = PredictedProfile,
        occupancyProfileParameters = occupancyProfileParameters,method="all")


    #Extracting Optimal matrix from Profile AccuracyEstimate

    opti<-.optimalExtraction(genomicProfileParameters,occupancyProfileParameters,ProfileAccuracy,optimalMethod)

    #return(list(opti,Occupancy,PredictedProfile,ProfileAccuracy))
    return(opti)
}
