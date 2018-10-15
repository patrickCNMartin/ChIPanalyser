##############################################
#############      S4 Methods    #############
##############################################

### Genomic Profile Parameters Accesor and setter methods
setMethod("PositionWeightMatrix","genomicProfileParameters",
    function(object) object@PWM
)


setReplaceMethod("PositionWeightMatrix",
    signature(object="genomicProfileParameters",
        value="matrix"),
    function(object, value){
        object@PWM<-value
    return(object)
    }
)


setMethod("PositionFrequencyMatrix","genomicProfileParameters",
    function(object) object@PFM
)

setReplaceMethod("PositionFrequencyMatrix",
    signature(object="genomicProfileParameters",
        value="matrix"),
    function(object, value){
        object@PFM<-value
    return(object)
    }
)
setReplaceMethod("PositionFrequencyMatrix",
    signature(object="genomicProfileParameters",
        value="character"),
    function(object,value){
        value <- .parsePFM(value,object@PFMFormat)
        object@PFM<-value
    return(object)
    }
)

setMethod("PFMFormat","genomicProfileParameters",
    function(object) object@PFMFormat
)

setReplaceMethod("PFMFormat",
    signature(object="genomicProfileParameters",
        value="character"),
    function(object, value){
        object@PFMFormat<-value
    return(object)
    }
)

setMethod("ScalingFactorPWM", "genomicProfileParameters",
    function(object) object@ScalingFactorPWM
)

setReplaceMethod("ScalingFactorPWM",
    signature(object="genomicProfileParameters",
        value="vector"),
    function(object,value){
        object@ScalingFactorPWM<-value
    return(object)
    }
)

setMethod("noOfSites", "genomicProfileParameters",
    function(object) object@noOfSites
)

setReplaceMethod("noOfSites",
    signature(object="genomicProfileParameters",
        value="vector"),
    function(object,value){
        object@noOfSites<-value
    return(object)
    }
)

setMethod("PWMpseudocount","genomicProfileParameters",
    function(object) object@PWMpseudocount
)

setReplaceMethod("PWMpseudocount",
    signature(object="genomicProfileParameters",
        value="numeric"),
    function(object, value){
        object@PWMpseudocount<-value
    return(object)
    }
)

setMethod("BPFrequency","genomicProfileParameters",
    function(object) object@BPFrequency
)

setReplaceMethod("BPFrequency",
    signature(object="genomicProfileParameters",
        value="vector"),
    function(object, value){
        object@BPFrequency<-value
    return(object)
    }
)
setReplaceMethod("BPFrequency",
    signature(object="genomicProfileParameters",
        value="DNAStringSet"),
    function(object, value){
        value <- .computeBPFrequency(value)
        object@BPFrequency<-value
    return(object)
    }
)


setMethod("naturalLog","genomicProfileParameters",
    function(object) object@naturalLog
)

setReplaceMethod("naturalLog",
    signature(object="genomicProfileParameters",
        value="logical"),
    function(object, value){
        object@naturalLog<-value
    return(object)
    }
)

setMethod("minPWMScore","genomicProfileParameters",
    function(object) object@minPWMScore
)


setMethod("maxPWMScore","genomicProfileParameters",
    function(object) object@maxPWMScore)


setMethod("PWMThreshold","genomicProfileParameters",
    function(object) object@PWMThreshold
)

setReplaceMethod("PWMThreshold",
    signature(object="genomicProfileParameters",
        value="numeric"),
    function(object, value){
        object@PWMThreshold<-value
    return(object)
    }
)


setMethod("AllSitesAboveThreshold","genomicProfileParameters",
    function(object) object@AllSitesAboveThreshold
)



setMethod("DNASequenceLength","genomicProfileParameters",
    function(object) object@DNASequenceLength
)

setReplaceMethod("DNASequenceLength",
    signature(object="genomicProfileParameters",
        value="vector"),
    function(object, value){
        object@DNASequenceLength<-value
    return(object)
    }
)

setMethod("averageExpPWMScore","genomicProfileParameters",
    function(object) object@averageExpPWMScore
)


setMethod("strandRule","genomicProfileParameters",
    function(object) object@strandRule
)

setReplaceMethod("strandRule",
    signature(object="genomicProfileParameters",
        value="character"),
    function(object, value){
        object@strandRule<-value
    return(object)
    }
)


setMethod("whichstrand","genomicProfileParameters",
    function(object) object@whichstrand
)

setReplaceMethod("whichstrand",
    signature(object="genomicProfileParameters",
        value="character"),
    function(object, value){
        object@whichstrand<-value
    return(object)
    }
)


setMethod("NoAccess","genomicProfileParameters",
    function(object) object@NoAccess
)





## Occupancy Profile Paramters Accessor and setter Methods

setMethod("ploidy","occupancyProfileParameters",
    function(object) object@ploidy
)

setReplaceMethod("ploidy",
    signature(object="occupancyProfileParameters",
        value="numeric"),
    function(object, value){
        object@ploidy<-value
    return(object)
    }
)

setMethod("boundMolecules","occupancyProfileParameters",
    function(object) object@boundMolecules
)

setReplaceMethod("boundMolecules",
    signature(object="occupancyProfileParameters",
        value="vector"),
    function(object, value){
        object@boundMolecules<-value
    return(object)
    }
)

setMethod("maxSignal","occupancyProfileParameters",
    function(object) object@maxSignal
)

setReplaceMethod("maxSignal",
    signature(object="occupancyProfileParameters",
        value="numeric"),
    function(object, value){
        object@maxSignal<-value
    return(object)
    }
)

setMethod("backgroundSignal","occupancyProfileParameters",
    function(object) object@backgroundSignal
)

setReplaceMethod("backgroundSignal",
    signature(object="occupancyProfileParameters",
        value="numeric"),
    function(object, value){
        object@backgroundSignal<-value
    return(object)
    }
)

setMethod("chipMean","occupancyProfileParameters",
    function(object) object@chipMean
)

setReplaceMethod("chipMean",
    signature(object="occupancyProfileParameters",
        value="numeric"),
    function(object, value){
        object@chipMean<-value
    return(object)
    }
)

setMethod("chipSd","occupancyProfileParameters",
    function(object) object@chipSd
)

setReplaceMethod("chipSd",
    signature(object="occupancyProfileParameters",
        value="numeric"),
    function(object, value){
        object@chipSd<-value
    return(object)
    }
)

setMethod("chipSmooth","occupancyProfileParameters",
    function(object) object@chipSmooth
)

setReplaceMethod("chipSmooth",
    signature(object="occupancyProfileParameters",
        value="vector"),
    function(object, value){
        object@chipSmooth<-value
    return(object)
    }
)

setMethod("removeBackground","occupancyProfileParameters",
    function(object) object@removeBackground
)

setReplaceMethod("removeBackground",
    signature(object="occupancyProfileParameters",
        value="vector"),
    function(object, value){
        object@removeBackground<-value
    return(object)
    }
)



setMethod("stepSize","occupancyProfileParameters",
    function(object) object@stepSize
)

setReplaceMethod("stepSize",
    signature(object="occupancyProfileParameters",
        value="numeric"),
    function(object, value){
        object@stepSize<-value
    return(object)
    }
)
