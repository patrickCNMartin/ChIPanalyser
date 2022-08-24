##############################################
#############      S4 Methods    #############
##############################################

### Genomic Profiles internal MEthods
setMethod("PositionWeightMatrix","genomicProfilesInternal",
    function(object) object@PWM
)


setReplaceMethod("PositionWeightMatrix",
    signature(object="genomicProfilesInternal",
        value="matrix"),
    function(object, value){
        object@PWM<-value
    return(object)
    }
)


setMethod("PositionFrequencyMatrix","genomicProfilesInternal",
    function(object) object@PFM
)

setReplaceMethod("PositionFrequencyMatrix",
    signature(object="genomicProfilesInternal",
        value="matrix"),
    function(object, value){
        object@PFM<-value
    return(object)
    }
)


setReplaceMethod("PositionFrequencyMatrix",
    signature(object="genomicProfilesInternal",
        value="character"),
    function(object,value){
        value <- .parsePFM(value,object@PFMFormat)
        object@PFM<-value
    return(object)
    }
)

setMethod("PFMFormat","genomicProfilesInternal",
    function(object) object@PFMFormat
)

setReplaceMethod("PFMFormat",
    signature(object="genomicProfilesInternal",
        value="character"),
    function(object, value){
        object@PFMFormat<-value
    return(object)
    }
)




setMethod("BPFrequency","genomicProfilesInternal",
    function(object) object@BPFrequency
)

setReplaceMethod("BPFrequency",
    signature(object="genomicProfilesInternal",
        value="vector"),
    function(object, value){
        object@BPFrequency<-value
    return(object)
    }
)


setReplaceMethod("BPFrequency",
    signature(object="genomicProfilesInternal",
        value="DNAStringSet"),
    function(object, value){
        value <- .computeBPFrequency(value)
        object@BPFrequency<-value
    return(object)
    }
)

setMethod("naturalLog","parameterOptions",
    function(object) object@naturalLog
)

setReplaceMethod("naturalLog",
    signature(object="parameterOptions",
        value="logical"),
    function(object, value){
        object@naturalLog<-value
    return(object)
    }
)


setMethod("minPWMScore","genomicProfilesInternal",
    function(object) object@minPWMScore
)
setReplaceMethod(".minPWMScore", signature(object="genomicProfilesInternal",
       value="vector"),
       function(object,value){
           object@minPWMScore<-value
           return(object)
         }
  )

setMethod("maxPWMScore","genomicProfilesInternal",
    function(object) object@maxPWMScore)

setReplaceMethod(".maxPWMScore", signature(object="genomicProfilesInternal",
    value="vector"),
       function(object,value){

           object@maxPWMScore<-value
           return(object)
         }
  )
setMethod("profiles","genomicProfilesInternal",
    function(object) object@profiles
)
setReplaceMethod(".profiles", signature(object="genomicProfilesInternal",
     value="GRList"),
     function(object,value){
       object@profiles<-value
       return(object)
     }
)


setMethod("DNASequenceLength","genomicProfilesInternal",
    function(object) object@DNASequenceLength
)

setReplaceMethod(".DNASequenceLength",
    signature(object="genomicProfilesInternal",
        value="vector"),
    function(object, value){
        object@DNASequenceLength<-value
    return(object)
    }
)

setMethod("averageExpPWMScore","genomicProfilesInternal",
    function(object) object@averageExpPWMScore
)

setReplaceMethod(".averageExpPWMScore",signature(object= "genomicProfilesInternal",
     value="numeric"),
    function(object,value){
       object@averageExpPWMScore<-value
       return(object)
    }
)

setMethod("drop","genomicProfilesInternal",
    function(object) object@drop
)

setReplaceMethod(".drop",signature(object= "genomicProfilesInternal",
    value="vector"),
    function(object, value){
       object@drop <- value
       return(object)
    }
)

setMethod(".tags","genomicProfilesInternal",
    function(object) object@tags
)

setReplaceMethod(".tags",
    signature(object="genomicProfilesInternal",
        value="character"),
    function(object, value){

        object@tags<-value

    return(object)
    }
)




## parameterOption methods

setMethod("PWMThreshold","parameterOptions",
    function(object) object@PWMThreshold
)

setReplaceMethod("PWMThreshold",
    signature(object="parameterOptions",
        value="numeric"),
    function(object, value){
        object@PWMThreshold<-value
    return(object)
    }
)


setMethod("strandRule","parameterOptions",
    function(object) object@strandRule
)

setReplaceMethod("strandRule",
    signature(object="parameterOptions",
        value="character"),
    function(object, value){
        object@strandRule<-value
    return(object)
    }
)


setMethod("whichstrand","parameterOptions",
    function(object) object@whichstrand
)

setReplaceMethod("whichstrand",
    signature(object="parameterOptions",
        value="character"),
    function(object, value){
        object@whichstrand<-value
    return(object)
    }
)




setMethod("lambdaPWM", "parameterOptions",
    function(object) object@lambdaPWM
)

setReplaceMethod("lambdaPWM",
    signature(object="parameterOptions",
        value="vector"),
    function(object,value){
        object@lambdaPWM<-value
    return(object)
    }
)

setMethod("noOfSites", "parameterOptions",
    function(object) object@noOfSites
)

setReplaceMethod("noOfSites",
    signature(object="parameterOptions",
        value="numeric"),
    function(object,value){
        object@noOfSites<-value
    return(object)
    }
)
setReplaceMethod("noOfSites",
    signature(object="parameterOptions",
        value="character"),
    function(object,value){
        object@noOfSites<-value
    return(object)
    }
)

setMethod("PWMpseudocount","parameterOptions",
    function(object) object@PWMpseudocount
)

setReplaceMethod("PWMpseudocount",
    signature(object="parameterOptions",
        value="numeric"),
    function(object, value){
        object@PWMpseudocount<-value
    return(object)
    }
)


setMethod("ploidy","parameterOptions",
    function(object) object@ploidy
)

setReplaceMethod("ploidy",
    signature(object="parameterOptions",
        value="numeric"),
    function(object, value){
        object@ploidy<-value
    return(object)
    }
)

setMethod("boundMolecules","parameterOptions",
    function(object) object@boundMolecules
)

setReplaceMethod("boundMolecules",
    signature(object="parameterOptions",
        value="vector"),
    function(object, value){
        object@boundMolecules<-value
    return(object)
    }
)

setMethod("maxSignal","parameterOptions",
    function(object) object@maxSignal
)

setReplaceMethod("maxSignal",
    signature(object="parameterOptions",
        value="numeric"),
    function(object, value){
        object@maxSignal<-value
    return(object)
    }
)

setMethod("backgroundSignal","parameterOptions",
    function(object) object@backgroundSignal
)

setReplaceMethod("backgroundSignal",
    signature(object="parameterOptions",
        value="numeric"),
    function(object, value){
        object@backgroundSignal<-value
    return(object)
    }
)

setMethod("chipMean","parameterOptions",
    function(object) object@chipMean
)

setReplaceMethod("chipMean",
    signature(object="parameterOptions",
        value="numeric"),
    function(object, value){
        object@chipMean<-value
    return(object)
    }
)

setMethod("chipSd","parameterOptions",
    function(object) object@chipSd
)

setReplaceMethod("chipSd",
    signature(object="parameterOptions",
        value="numeric"),
    function(object, value){
        object@chipSd<-value
    return(object)
    }
)

setMethod("chipSmooth","parameterOptions",
    function(object) object@chipSmooth
)

setReplaceMethod("chipSmooth",
    signature(object="parameterOptions",
        value="vector"),
    function(object, value){
        object@chipSmooth<-value
    return(object)
    }
)

setMethod("removeBackground","parameterOptions",
    function(object) object@removeBackground
)

setReplaceMethod("removeBackground",
    signature(object="parameterOptions",
        value="vector"),
    function(object, value){
        object@removeBackground<-value
    return(object)
    }
)


setMethod(".ZeroBackground","parameterOptions",
    function(object) object@ZeroBackground
)

setReplaceMethod(".ZeroBackground",
    signature(object="parameterOptions",
        value="vector"),
    function(object, value){
        object@ZeroBackground<-value
    return(object)
    }
)


setMethod("stepSize","parameterOptions",
    function(object) object@stepSize
)

setReplaceMethod("stepSize",
    signature(object="parameterOptions",
        value="numeric"),
    function(object, value){
        object@stepSize<-value
    return(object)
    }
)

setMethod("lociWidth","parameterOptions",
    function(object) object@lociWidth
)



setReplaceMethod("lociWidth",
    signature(object="parameterOptions",
        value="numeric"),
    function(object, value){
        object@lociWidth<-value
    return(object)
    }
)
setReplaceMethod("noiseFilter",
    signature(object="parameterOptions",
        value="character"),
    function(object, value){

        object@noiseFilter<-value

    return(object)
    }
)

setMethod("noiseFilter","parameterOptions",
    function(object) object@noiseFilter
)



setMethod(".paramTag","parameterOptions",
    function(object) object@paramTag
)

setReplaceMethod(".paramTag",
    signature(object="parameterOptions",
        value="character"),
    function(object, value){

        object@paramTag<-value

         return(object)
    }
)


setMethod("loci","ChIPScore",
    function(object) object@loci
)

setReplaceMethod(".loci",
    signature(object="ChIPScore",
        value="loci"),
    function(object, value){

        object@loci<-value

         return(object)
    }
)

setMethod("scores","ChIPScore",
    function(object) object@scores
)

setReplaceMethod(".scores",
    signature(object="ChIPScore",
        value="list"),
    function(object, value){

        object@scores<-value

         return(object)
    }
)
