##########################################
#############  S4 Generics   #############
##########################################



#### Genomic Profile Parameters ####

#setGeneric(".initialize",
    #function(object) standardGeneric(".initialize"))

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

setGeneric("lambdaPWM",
    function(object) standardGeneric("lambdaPWM"))

setGeneric("lambdaPWM<-",
    function(object, value) standardGeneric("lambdaPWM<-"))

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
setGeneric(".minPWMScore<-",
      function(object,value)standardGeneric(".minPWMScore<-"))

setGeneric("maxPWMScore",
    function(object) standardGeneric("maxPWMScore"))
setGeneric(".maxPWMScore<-",
     function(object,value)standardGeneric(".maxPWMScore<-"))

setGeneric(".tags",
      function(object) standardGeneric(".tags"))
setGeneric(".tags<-",
      function(object,value) standardGeneric(".tags<-"))

setGeneric("PWMThreshold",
    function(object) standardGeneric("PWMThreshold"))

setGeneric("PWMThreshold<-",
    function(object, value) standardGeneric("PWMThreshold<-"))

setGeneric("profiles",
    function(object) standardGeneric("profiles"))
setGeneric(".profiles<-",function(object,value)standardGeneric(".profiles<-"))



setGeneric("DNASequenceLength",
    function(object) standardGeneric("DNASequenceLength"))

setGeneric(".DNASequenceLength<-",
    function(object, value) standardGeneric(".DNASequenceLength<-"))

setGeneric("averageExpPWMScore",
    function(object) standardGeneric("averageExpPWMScore"))
setGeneric(".averageExpPWMScore<-",
    function(object,value) standardGeneric(".averageExpPWMScore<-"))


setGeneric("strandRule",
    function(object) standardGeneric("strandRule"))

setGeneric("strandRule<-",
    function(object,value) standardGeneric("strandRule<-"))

setGeneric("whichstrand",
    function(object) standardGeneric("whichstrand"))

setGeneric("whichstrand<-",
    function(object,value) standardGeneric("whichstrand<-"))



setGeneric(".ZeroBackground",
    function(object) standardGeneric(".ZeroBackground"))
setGeneric(".ZeroBackground<-",
    function(object,value) standardGeneric(".ZeroBackground<-"))

setGeneric(".drop<-",
    function(object,value) standardGeneric(".drop<-"))

setGeneric("drop",
    function(object) standardGeneric("drop"))




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

setGeneric("lociWidth",
        function(object) standardGeneric("lociWidth"))

setGeneric("lociWidth<-",
        function(object, value) standardGeneric("lociWidth<-"))

setGeneric("noiseFilter",
          function(object) standardGeneric("noiseFilter"))

setGeneric("noiseFilter<-",
        function(object, value) standardGeneric("noiseFilter<-"))

setGeneric(".paramTag",
      function(object) standardGeneric(".paramTag"))
setGeneric(".paramTag<-",
      function(object,value) standardGeneric(".paramTag<-"))


setGeneric("loci",
      function(object) standardGeneric("loci"))
setGeneric("loci<-",
      function(object,value) standardGeneric("loci<-"))

setGeneric("scores",
      function(object) standardGeneric("scores"))
setGeneric("scores<-",
      function(object,value) standardGeneric("scores<-"))
