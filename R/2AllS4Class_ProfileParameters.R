###############################################################################
#########################     All s4 class ####################################
###############################################################################

## Class Union for Genomic Profile internal
setClassUnion("GRList",c("CompressedGRangesList","list","GRanges"))
setClassUnion("nos",c("numeric","character"))

# apparently PFM does not come up so we will try this
# this works dont know why though 
# do not touch this
utils::globalVariables("PFM")

## Genomic Profile Parameters
setClass("genomicProfilesInternal",
    slots = c(PWM = "matrix",
    PFM = "matrix",
    PFMFormat = "character",
    BPFrequency = "vector",
    minPWMScore = "vector",
    maxPWMScore = "vector",
    profiles = "GRList",
    DNASequenceLength = "vector",
    averageExpPWMScore = "vector",
    ZeroBackground = "vector",
    drop="vector",
    tags="character"),

    prototype = prototype(PFMFormat="raw",
    BPFrequency = rep(0.25,4),
    tags="empty"),

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



    return(TRUE)
    }
)


#setClass("ChIPScore",contains="VIRTUAL")
setClass("parameterOptions",

    slots = c(ploidy = "numeric",
    boundMolecules = "vector",
    backgroundSignal ="numeric",
    maxSignal = "numeric",
    lociWidth ="numeric",
    chipMean = "numeric",
    chipSd = "numeric",
    chipSmooth = "vector",
    stepSize = "numeric",
    removeBackground = "numeric",
    noiseFilter="character",
    PWMThreshold="numeric",
    strandRule="character",
    whichstrand="character",
    lambdaPWM="vector",
    naturalLog="logical",
    noOfSites="nos",
    PWMpseudocount="numeric",
    paramTag="character"),

    prototype = prototype(ploidy = 2 ,
    boundMolecules = 1000 ,
    backgroundSignal = 0 ,
    maxSignal = 1 ,
    lociWidth = 20000,
    chipMean = 200 ,
    chipSd = 200 ,
    chipSmooth = 250 ,
    stepSize = 10 ,
    removeBackground = 0,
    noiseFilter="zero",
    naturalLog=TRUE,
    PWMThreshold=0.7,
    strandRule="max",
    whichstrand="+-",
    PWMpseudocount=1,
    lambdaPWM=1,
    noOfSites="all",
    paramTag="empty"),

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
    if(object@lociWidth < 0){
    stop("Occupancy Profile Parameters: lociWidth  must be positive.")
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
    if(class(object@removeBackground)!="numeric"){
    stop("removeBackground Must be numeric")
    }
    if(any(object@lambdaPWM < 0)){
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
    if(all(object@noiseFilter %in% c("zero","mean","median","sigmoid")) == FALSE){
    stop("Genomic Profile Paramters:
    noiseFilter sould be one of the following:
    zero, mean,median, or sigmoid")
    }
    return(TRUE)
    }
)

setClassUnion("loci",c("GRanges","NULL"))
setClass("ChIPScore",
     contains="parameterOptions",
    slots=c(scores="list",
            loci="loci" ),
       validity=function(object){

        if(class(object@scores)!="list"){
           stop("Internal Error: ChIPScore score is not a list")
        }

        if(!is.null(object@loci)& class(object@loci)!="GRanges"){
           stop("Internal Error: loci is not empty but it ain't no GRanges")
        }

        return(TRUE)
    }
)


setClass("genomicProfiles",
         contains=c("genomicProfilesInternal","parameterOptions")
)



.genomicProfilesInternalContructor <- function(
    PWM=NULL ,
    PFM=NULL ,
    PFMFormat="raw",
    BPFrequency = rep(0.25,4)
    ){
    if(!is.null(PWM)| !is.null(PFM)){
       tags <- "PWM_stage"
    } else{
       tags <- "empty"
    }
    return(new("genomicProfilesInternal",
        PWM = PWM,
        PFM = PFM,
        PFMFormat =PFMFormat,
        BPFrequency = BPFrequency,
        tags = tags
        ))


}

genomicProfiles <- function(..., parameterOptions=NULL,genomicProfiles=NULL,ChIPScore=NULL){
      args<-list(...)


    ## seperating argumants for each S4

    PO<-c("ploidy","boundMolecules","backgroundSignal","maxSignal","lociWidth",
              "chipMean","chipSd","chipSmooth","stepSize","removeBackground","noiseFilter",
               "naturalLog","noOfSites","PWMThreshold","strandRule","whichstrand",
               "PWMpseudocount","lambdaPWM")
    GP<-c("PWM","PFM","PFMFormat","BPFrequency")
    POargs<-args[names(args) %in% PO]
    GPargs<-args[names(args) %in% GP]
    ## creating argument parsing vector
    GPLocal<-match(GP,names(GPargs))
    POLocal<-match(PO,names(POargs))
    localEnvir<-environment()
    ## assigning values to variable for argument parsing

    mapply(function(GPLocal,GP,localEnvir,GPargs){
         if(is.na(GPLocal)){
         assign(GP,NULL,envir=localEnvir)
       }else{
         assign(GP,GPargs[[GPLocal]],envir=localEnvir)
       }},GPLocal,GP,MoreArgs=list(localEnvir,GPargs))

    mapply(function(POLocal,PO,localEnvir,POargs){
        if(is.na(POLocal)){
          assign(PO,NULL,envir=localEnvir)
        }else{
          assign(PO,POargs[[POLocal]],envir=localEnvir)
        }},POLocal,PO,MoreArgs=list(localEnvir,POargs))
     # tag change
     if(!is.null(PWM)| !is.null(PFM)){
         tags <- "PWM_stage"
      } else{
         tags <- "empty"
      }
       return(new("genomicProfiles",ploidy = ploidy ,
       boundMolecules = boundMolecules ,
       backgroundSignal = backgroundSignal ,
       maxSignal = maxSignal ,
       lociWidth = lociWidth,
       chipMean = chipMean ,
       chipSd = chipSd ,
       chipSmooth = chipSmooth ,
       stepSize = stepSize ,
       removeBackground = removeBackground,
       noiseFilter = noiseFilter,
       naturalLog = naturalLog,
       noOfSites = noOfSites,
       PWMThreshold = PWMThreshold,
       strandRule = strandRule,
       whichstrand = whichstrand,
       PWMpseudocount = PWMpseudocount,
       lambdaPWM = lambdaPWM,
       PWM = PWM,
       PFM = PFM,
       PFMFormat =PFMFormat,
       BPFrequency = BPFrequency,
       tags= tags))
}




parameterOptions <- function(
  ploidy = 2 ,
  boundMolecules = 1000 ,
  backgroundSignal = 0 ,
  maxSignal = 1 ,
  lociWidth = 20000,
  chipMean = 200 ,
  chipSd = 200 ,
  chipSmooth = 250 ,
  stepSize = 10 ,
  removeBackground = 0,
  noiseFilter="zero",
  naturalLog=TRUE,
  noOfSites = "all",
  PWMThreshold=0.7,
  strandRule="max",
  whichstrand="+-",
  PWMpseudocount=1,
  lambdaPWM=1){

     ## step check


    PO<-new("parameterOptions",
        ploidy = ploidy ,
        boundMolecules = boundMolecules ,
        backgroundSignal = backgroundSignal ,
        maxSignal = maxSignal ,
        lociWidth = lociWidth,
        chipMean = chipMean ,
        chipSd = chipSd ,
        chipSmooth = chipSmooth ,
        stepSize = stepSize ,
        removeBackground = removeBackground,
        noiseFilter = noiseFilter,
        naturalLog = naturalLog,
        noOfSites = noOfSites,
        PWMThreshold = PWMThreshold,
        strandRule = strandRule,
        whichstrand = whichstrand,
        PWMpseudocount = PWMpseudocount,
        lambdaPWM = lambdaPWM
        )
    return(PO)
}

.ChIPScore <-function(scores=NULL,loci=NULL,maxSignal=NULL,
                      backgroundSignal=NULL,lociWidth=NULL,paramTag="empty"){
    return(new("ChIPScore",
        scores=scores,
        loci=loci,
        maxSignal=maxSignal,
        backgroundSignal=backgroundSignal,
        lociWidth=lociWidth,
        paramTag=paramTag))
}
