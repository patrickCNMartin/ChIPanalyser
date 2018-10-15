##########################################
#############  S4 Generics   #############
##########################################



#### Genomic Profile Parameters ####

setGeneric(".initialize",
    function(object) standardGeneric(".initialize"))

setGeneric(".generatePWM",
    function(object) standardGeneric(".generatePWM"))

setGeneric("PositionWeightMatrix",
    function(object) standardGeneric("PositionWeightMatrix"))

setGeneric("PositionWeightMatrix<-",
    function(object, value) standardGeneric("PositionWeightMatrix<-"))

setGeneric("PositionFrequencyMatrix",
    function(object) standardGeneric("PositionFrequencyMatrix"))

setGeneric("PositionFrequencyMatrix<-",
    function(object, value) standardGeneric("PositionFrequencyMatrix<-"))

setGeneric("PFMFormat",
    function(object) standardGeneric("PFMFormat"))

setGeneric("PFMFormat<-",
    function(object, value) standardGeneric("PFMFormat<-"))

setGeneric("ScalingFactorPWM",
    function(object) standardGeneric("ScalingFactorPWM"))

setGeneric("ScalingFactorPWM<-",
    function(object, value) standardGeneric("ScalingFactorPWM<-"))

setGeneric("noOfSites",
    function(object) standardGeneric("noOfSites"))

setGeneric("noOfSites<-",
    function(object, value) standardGeneric("noOfSites<-"))

setGeneric("PWMpseudocount",
    function(object) standardGeneric("PWMpseudocount"))

setGeneric("PWMpseudocount<-",
    function(object, value) standardGeneric("PWMpseudocount<-"))

setGeneric("BPFrequency",
    function(object) standardGeneric("BPFrequency"))

setGeneric("BPFrequency<-",
    function(object, value) standardGeneric("BPFrequency<-"))

setGeneric("naturalLog",
    function(object) standardGeneric("naturalLog"))

setGeneric("naturalLog<-",
    function(object, value) standardGeneric("naturalLog<-"))

setGeneric("minPWMScore",
    function(object) standardGeneric("minPWMScore"))

setGeneric("maxPWMScore",
    function(object) standardGeneric("maxPWMScore"))



setGeneric("PWMThreshold",
    function(object) standardGeneric("PWMThreshold"))

setGeneric("PWMThreshold<-",
    function(object, value) standardGeneric("PWMThreshold<-"))

setGeneric("AllSitesAboveThreshold",
    function(object) standardGeneric("AllSitesAboveThreshold"))



setGeneric("DNASequenceLength",
    function(object) standardGeneric("DNASequenceLength"))

setGeneric("DNASequenceLength<-",
    function(object, value) standardGeneric("DNASequenceLength<-"))

setGeneric("averageExpPWMScore",
    function(object) standardGeneric("averageExpPWMScore"))



setGeneric("strandRule",
    function(object) standardGeneric("strandRule"))

setGeneric("strandRule<-",
    function(object,value) standardGeneric("strandRule<-"))

setGeneric("whichstrand",
    function(object) standardGeneric("whichstrand"))

setGeneric("whichstrand<-",
    function(object,value) standardGeneric("whichstrand<-"))

setGeneric("NoAccess",
    function(object) standardGeneric("NoAccess"))




#### Occupancy Profile Parameters ####




setGeneric("ploidy",
    function(object) standardGeneric("ploidy"))

setGeneric("ploidy<-",
    function(object, value) standardGeneric("ploidy<-"))

setGeneric("boundMolecules",
    function(object) standardGeneric("boundMolecules"))

setGeneric("boundMolecules<-",
    function(object, value) standardGeneric("boundMolecules<-"))

setGeneric("maxSignal",
    function(object) standardGeneric("maxSignal"))

setGeneric("maxSignal<-",
    function(object, value) standardGeneric("maxSignal<-"))

setGeneric("backgroundSignal",
    function(object) standardGeneric("backgroundSignal"))

setGeneric("backgroundSignal<-",
    function(object, value) standardGeneric("backgroundSignal<-"))

setGeneric("chipMean",
    function(object) standardGeneric("chipMean"))

setGeneric("chipMean<-",
    function(object, value) standardGeneric("chipMean<-"))

setGeneric("chipSd",
    function(object) standardGeneric("chipSd"))

setGeneric("chipSd<-",
    function(object, value) standardGeneric("chipSd<-"))

setGeneric("chipSmooth",
    function(object) standardGeneric("chipSmooth"))

setGeneric("chipSmooth<-",
    function(object, value) standardGeneric("chipSmooth<-"))

setGeneric("removeBackground",
    function(object) standardGeneric("removeBackground"))

setGeneric("removeBackground<-",
    function(object, value) standardGeneric("removeBackground<-"))


setGeneric("stepSize",
    function(object) standardGeneric("stepSize"))

setGeneric("stepSize<-",
    function(object, value) standardGeneric("stepSize<-"))
