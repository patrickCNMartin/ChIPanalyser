#############################################
###############  S4 Objects   ###############
#############################################

## Occupancy Profile Parameters

setClass("occupancyProfileParameters",
    slots = c(ploidy = "numeric",
    boundMolecules = "vector",
    backgroundSignal ="numeric",
    maxSignal = "numeric",
    chipMean = "numeric",
    chipSd = "numeric",
    chipSmooth = "vector",
    stepSize = "numeric",
    removeBackground = "numeric",
    thetaThreshold = "numeric"),

    prototype = prototype(ploidy = 2 ,
    boundMolecules = 2000 ,
    backgroundSignal = 0 ,
    maxSignal = 1 ,
    chipMean = 200 ,
    chipSd = 200 ,
    chipSmooth = 250 ,
    stepSize = 10 ,
    removeBackground = 0 ,
    thetaThreshold = 0.1 ),

    validity = function(object){
    if(object@ploidy < 0){
    stop("Occupancy Profile Parameters: Ploidy must be positive.")
    }
    if(any(object@boundMolecules < 0)){
    stop("Occupancy Profile Parameters:
    Bound Molecules must be positive.")
    }
    if(object@backgroundSignal < 0){
    stop("Occupancy Profile Parameters:
    backgroundSignal must be positive")
    }
    if(object@maxSignal < 0 |object@maxSignal < object@backgroundSignal ){
    stop("Occupancy Profile Parameters:
    maxSignal must be positive and greater then backgroundSignal")
    }
    if(object@chipMean < 0){
    stop("Occupancy Profile Parameters: ChIP Mean must be positive.")
    }
    if(object@chipSd < 0){
    stop("Occupancy Profile Parameters: ChIP SD must be positive.")
    }
    if(object@chipSmooth <= 0){
    stop("Occupancy Profile Parameters: ChIP chipSmooth must be positive
    or equal to zero.")
    }
    if(object@stepSize < 0){
    stop("Occupancy Profile Parameters: Step Size must be positive.")
    }
    if(object@removeBackground < 0){
    stop("Occupancy Profile Parameters:
    removeBackground must be positive")
    }
    if(object@thetaThreshold < 0){
    stop("Occupancy Profile Parameters:
    thetaThreshold must be positive")
    }
    return(TRUE)
    }
)

## Class Union for Genomic Profile Parameters slot
setClassUnion("GRList",c("GRangesList","list"))

## Genomic Profile Parameters
setClass("genomicProfileParameters",
    slots = c(PWM = "matrix",
    PFM = "matrix",
    PFMFormat = "character",
    ScalingFactorPWM = "vector",
    PWMpseudocount = "numeric",
    BPFrequency = "vector",
    naturalLog = "logical",
    noOfSites = "vector",
    minPWMScore = "vector",
    maxPWMScore = "vector",
    PWMThreshold = "numeric",
    AllSitesAboveThreshold = "GRList",
    DNASequenceLength = "vector",
    averageExpPWMScore = "vector",
    strandRule = "character",
    whichstrand = "character",
    NoAccess = "vector",
    ZeroBackground = "vector"),

    prototype = prototype(PFMFormat="raw",
    ScalingFactorPWM = 1.25,
    PWMpseudocount = 1,
    noOfSites = 0,
    BPFrequency = rep(0.25,4),
    naturalLog = FALSE,
    PWMThreshold = 0.7,
    strandRule = "max",
    whichstrand = "+-"),

    contains = "occupancyProfileParameters",

    validity= function(object){
    if(length(object@PFM) < 1 &
        length(object@PWM) > 0 &
        all(object@PWM >= 0)){
    stop("Genomic Profile Parameters Validation:
    PWM does not seem to be a Position Weight Matrix.")
    }


    if(all(object@PFMFormat %in% c("raw","transfac","JASPAR","sequences","matrix"))
    ==FALSE ){
        stop("Genomic Profile Parameters Validation:
        PFM file is not one of the following formats:
        raw,transfac,JASPAR,sequences or matrix.")
    }
    if(any(object@ScalingFactorPWM < 0)){
    stop("Genomic Profile Parameters Validation:
    Lambda must be positive.")
    }
    if(object@PWMpseudocount < 0 ){
    stop("Genomic Profile Parameters: psuedocount must be positive")
    }
    if(object@PWMThreshold < 0 | object@PWMThreshold >= 1){
    stop("Genomic Profile Parameters Validation:
    PWMThreshold must be between 0 and 1")
    }
    if(any(object@BPFrequency < 0)){
    stop("Genomic Profile Parameters:
    Base Pair frequency contains negative frequencies")
    }
    if(!is.logical(object@naturalLog)){
    stop("Genomic Profile Parameters: naturalLog should either be
    TRUE or FALSE")
    }
    if(all(object@strandRule %in% c("max","mean","sum")) == FALSE){
    stop("Genomic Profile Paramters:
    strandRule should be one of the following:
    max, mean, or sum")
    }
    if(all(object@whichstrand %in% c("+","-","+-","-+")) == FALSE){
    stop("Genomic Profile Paramters:
    whichstrand sould be one of the following:
    +, -, +-, or -+")
    }
    return(TRUE)
    }
)


genomicProfileParameters <- function(
    PWM=NULL ,
    PFM=NULL ,
    PFMFormat="raw",
    ScalingFactorPWM = 1,
    PWMpseudocount = 1,
    noOfSites = 0,
    BPFrequency = rep(0.25,4),
    naturalLog = FALSE,
    PWMThreshold = 0.7,
    strandRule = "max",
    whichstrand = "+-"){
    return(new("genomicProfileParameters",
        PWM = PWM,
        PFM = PFM,
        PFMFormat =PFMFormat,
        ScalingFactorPWM = ScalingFactorPWM,
        PWMpseudocount = PWMpseudocount,
        noOfSites = noOfSites,
        BPFrequency = BPFrequency,
        naturalLog = naturalLog,
        PWMThreshold = PWMThreshold,
        strandRule = strandRule,
        whichstrand = whichstrand))


}


occupancyProfileParameters <- function(
    ploidy = 2 ,
    boundMolecules = 1000 ,
    backgroundSignal = 0 ,
    maxSignal = 1 ,
    chipMean = 150 ,
    chipSd = 150 ,
    chipSmooth = 250 ,
    stepSize = 10 ,
    removeBackground = 0 ,
    thetaThreshold = 0.1){
    return(new("occupancyProfileParameters",
        ploidy = ploidy ,
        boundMolecules = boundMolecules ,
        backgroundSignal = backgroundSignal ,
        maxSignal = maxSignal ,
        chipMean = chipMean ,
        chipSd = chipSd ,
        chipSmooth = chipSmooth ,
        stepSize = stepSize ,
        removeBackground = removeBackground ,
        thetaThreshold = thetaThreshold))


}
