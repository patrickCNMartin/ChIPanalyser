##############################################
#############      S4 Methods    #############
##############################################



### initialize Genomic Profile by generating a PWM from a PFM
setMethod(".generatePWM","genomicProfilesInternal",
    function(object){
        if(object@noOfSites=="all"){
            object@noOfSites<-0
        }
        newPWM <- .computePWM(object@PFM,
            object@PWMpseudocount,
            object@BPFrequency,
            object@naturalLog,
            object@noOfSites)
        rownames(newPWM) <- c("A", "C", "G", "T")
        object@PWM <- newPWM
        return(object)
    }
)


setMethod(f="initialize",
    signature = "parameterOptions",
    definition = function(.Object,
        ploidy = .Object@ploidy,
        boundMolecules = .Object@boundMolecules,
        backgroundSignal = .Object@backgroundSignal,
        maxSignal = .Object@maxSignal,
        lociWidth = .Object@lociWidth,
        chipMean = .Object@chipMean,
        chipSd = .Object@chipSd,
        chipSmooth = .Object@chipSmooth,
        stepSize = .Object@stepSize,
        noiseFilter= .Object@noiseFilter,
        removeBackground = .Object@removeBackground,
        lambdaPWM=.Object@lambdaPWM,
        PWMpseudocount=.Object@PWMpseudocount,
        naturalLog=.Object@naturalLog,
        noOfSites=.Object@noOfSites,
        PWMThreshold=.Object@PWMThreshold,
        strandRule=.Object@strandRule,
        whichstrand=.Object@whichstrand
        ){

        if(length(ploidy)>0){
            if(class(ploidy)=="numeric"){
        NewPloidy<-ploidy
        .Object@ploidy <- NewPloidy
        validObject(.Object)
        } else {
        stop("ploidy is not numeric")}}
        if(length(lambdaPWM)>0){
            if(class(lambdaPWM)=="numeric"){
        NewLambda<-lambdaPWM
        .Object@lambdaPWM<-NewLambda
        validObject(.Object)
        } else {
        stop("lambdaPWM is not numeric
        or a vector of numeric values")
        }}
        if(length(boundMolecules)>0){
            if(class(boundMolecules)=="numeric"){
        NewboundMolecules<-boundMolecules
        .Object@boundMolecules <- NewboundMolecules
        validObject(.Object)
        } else {
        stop("boundMolecules is not numeric")}}

        if(length(backgroundSignal)>0){
            if(class(backgroundSignal)=="numeric"){
        NewbackgroundSignal<-backgroundSignal
        .Object@backgroundSignal <- NewbackgroundSignal
        validObject(.Object)
        } else {
        stop("backgroundSignal is not numeric")}}

        if(length(maxSignal)>0){
            if(class(maxSignal)=="numeric"){
        NewmaxSignal<-maxSignal
        .Object@maxSignal <- NewmaxSignal
        validObject(.Object)
        } else {
        stop("maxSignal is not numeric")}}

        if(length(lociWidth)>0){
            if(class(lociWidth)=="numeric"){
        NewmlociWidth<-lociWidth
        .Object@lociWidth <- NewmlociWidth
        validObject(.Object)
        } else {
        stop("maxSignal is not numeric")}}

        if(length(chipMean)>0){
            if(class(chipMean)=="numeric"){
        NewchipMean<-chipMean
        .Object@chipMean <- NewchipMean
        validObject(.Object)
        } else {
        stop("chipMean is not numeric")}}

        if(length(chipSd)>0){
            if(class(chipSd)=="numeric"){
        NewchipSd<-chipSd
        .Object@chipSd <- NewchipSd
        validObject(.Object)
        } else {
        stop("chipSd is not numeric")}}

        if(length(chipSmooth)>0){
            if(class(chipSmooth)=="numeric"){
        NewchipSmooth<-chipSmooth
        .Object@chipSmooth <- NewchipSmooth
        validObject(.Object)
        }else{
        stop("chipSmooth is not numeric")
        }}
        if(length(stepSize)>0){
            if(class(stepSize)=="numeric"){
        NewstepSize<-stepSize
        .Object@stepSize <- NewstepSize
        validObject(.Object)
        } else {
        stop("stepSize is not numeric")}}

        if(length(removeBackground)>0){
            if(class(removeBackground)=="numeric"){
        NewremoveBackground<-removeBackground
        .Object@removeBackground <- NewremoveBackground
        validObject(.Object)
        } else {
        stop("removeBackground is not numeric")}}
        if(length(PWMpseudocount)>0){
            if(class(PWMpseudocount)=="numeric"){
        Newpseudocount<-PWMpseudocount
        .Object@PWMpseudocount<-Newpseudocount
        validObject(.Object)
        } else {
        stop("PWMpseudocount is not numeric")
        }}

        if(length(naturalLog)>0 & is.logical(naturalLog)){
        NewnaturalLog<-naturalLog
        .Object@naturalLog<-NewnaturalLog
        validObject(.Object)
        }
        if(length(noOfSites)>0){
            if(class(noOfSites)=="numeric"|class(noOfSites)=="character"){
        NewnoOfSites<-noOfSites
        .Object@noOfSites<-NewnoOfSites
        validObject(.Object)
        } else {
        stop("noOfSites is not numeric or 'all' - for all sites")
        }}
        if(length(PWMThreshold)>0){
            if(class(PWMThreshold)=="numeric"){
        NewPWMThreshold<-PWMThreshold
        .Object@PWMThreshold<-NewPWMThreshold
        validObject(.Object)
        } else {
        stop("PWMThreshold is not numeric")
        }}

        if(length(strandRule)>0){
            if(class(strandRule)=="character"){
        NewStrandRule<-strandRule
        .Object@strandRule<-NewStrandRule
        validObject(.Object)
        } else {
        stop("strandRule is not a character string of one
        of the following: max, mean or sum")
        }}

        if(length(whichstrand)>0){
            if(class(whichstrand)=="character"){
        NewStrand<-whichstrand
        .Object@whichstrand<-NewStrand
        validObject(.Object)
        } else {
        stop("whichstrand is not a character string of one of the
        following: +-, -+ , - or +")
        }}
        if(length(noiseFilter)>0){
            if(class(noiseFilter)=="character"){

        .Object@noiseFilter<-noiseFilter
        validObject(.Object)
        } else {
        stop("noiseFilter is not a character string of one of the
        following: zero, mean, median or sigmoid")
        }}

    return(.Object)
    }
)



setMethod(f="initialize",
    signature = "ChIPScore",
    definition = function(.Object,
        scores = .Object@scores,
        loci = .Object@loci,
        maxSignal= .Object@maxSignal,
        backgroundSignal=.Object@backgroundSignal,
        lociWidth=.Object@lociWidth,
        paramTag=.Object@paramTag){

        if(length(scores)>0){

            .Object@scores <- scores

            validObject(.Object)
        }

        if(length(loci)>0){

            .Object@loci <- loci
            validObject(.Object)
        }
        if(length(maxSignal)>0){

            .Object@maxSignal <- maxSignal
            validObject(.Object)
        }
        if(length(backgroundSignal)>0){

            .Object@backgroundSignal <- backgroundSignal
            validObject(.Object)
        }
        if(length(paramTag)>0){

            .Object@paramTag <- paramTag
            validObject(.Object)
        }
        if(length(lociWidth)>0){

            .Object@lociWidth <- lociWidth
            validObject(.Object)
        }
        return(.Object)
    }
)

setMethod(f="initialize",
  signature = "genomicProfiles",
    definition = function(.Object,
    ploidy = .Object@ploidy,
    boundMolecules = .Object@boundMolecules,
    backgroundSignal = .Object@backgroundSignal,
    maxSignal = .Object@maxSignal,
    lociWidth = .Object@lociWidth,
    chipMean = .Object@chipMean,
    chipSd = .Object@chipSd,
    chipSmooth = .Object@chipSmooth,
    stepSize = .Object@stepSize,
    noiseFilter= .Object@noiseFilter,
    removeBackground = .Object@removeBackground,
    lambdaPWM=.Object@lambdaPWM,
    PWMpseudocount=.Object@PWMpseudocount,
    naturalLog=.Object@naturalLog,
    noOfSites=.Object@noOfSites,
    PWMThreshold=.Object@PWMThreshold,
    strandRule=.Object@strandRule,
    whichstrand=.Object@whichstrand,
    paramTags=.Object@paramTags,
    PFM=.Object@PFM,
    PWM=.Object@PWM,
    PFMFormat=.Object@PFMFormat,
    BPFrequency=.Object@BPFrequency,
    minPWMScore=.Object@minPWMScore,
    maxPWMScore=.Object@maxPWMScore,
    profiles = .Object@profiles,
    DNASequenceLength= .Object@DNASequenceLength,
    averageExpPWMScore=.Object@averageExpPWMScore,
    ZeroBackground=.Object@ZeroBackground,
    drop=.Object@drop,
    tags=.Object@tags){
      if(class(BPFrequency) != "vector" &
      (class(BPFrequency) == "BSgenome" |
      class(BPFrequency) == "DNAStringSet")){
      BPFrequency<-.computeBPFrequency(BPFrequency)
      .Object@BPFrequency<-BPFrequency
      } else if(!is.null(BPFrequency) & class(BPFrequency)=="vector") {
      NewBPFrequency<-BPFrequency
      .Object@BPFrequency<-NewBPFrequency
      validObject(.Object)
      }

      if(length(PFMFormat)>0){
          if(class(PFMFormat)=="character"){
      NewPFMFormat<-PFMFormat
      .Object@PFMFormat<-NewPFMFormat
      validObject(.Object)
      } else {
      stop("PFMFormat is not a character string of one of the following:
      raw, transfac, JASPAR, sequences or matrix")
      }}

      if(length(minPWMScore)>0){

      .Object@minPWMScore<-minPWMScore
      validObject(.Object)
      }

      if(length(maxPWMScore)>0){

      .Object@maxPWMScore<-maxPWMScore
      validObject(.Object)
      }
      if(!is.null(PFM)){
          if(class(PFM)!="matrix"){
      PFMmatrix<-.parsePFM(PFM,PFMFormat)
      .Object@PFM<-PFMmatrix
      validObject(.Object)
      } else {
      PFMMatrix<-PFM
      .Object@PFM<-PFMMatrix
      validObject(.Object)
      }}

      if(length(PFM)<1 & length(PWM)>0){
      NewPWM<-PositionWeightMatrix
      .Object@PWM<-NewPWM
      validObject(.Object)
      }
      if(length(PFM)>0 & length(PWM)<1){
      New<-.generatePWM(.Object)
      New<-New@PWM
      .Object@PWM<-New
      validObject(.Object)
      }
      if(length(profiles)>0){

      .Object@profiles<-profiles
      validObject(.Object)
      }

      if(length(DNASequenceLength)>0){

      .Object@DNASequenceLength<-DNASequenceLength
      validObject(.Object)
      }

      if(length(averageExpPWMScore)>0){

      .Object@averageExpPWMScore<-averageExpPWMScore
      validObject(.Object)
      }

      if(length(drop)>0){

      .Object@drop<-drop
      validObject(.Object)
      }

      if(length(tags)>0){

      .Object@tags<-tags
      validObject(.Object)
      }
      if(length(ploidy)>0){
          if(class(ploidy)=="numeric"){
      NewPloidy<-ploidy
      .Object@ploidy <- NewPloidy
      validObject(.Object)
      } else {
      stop("ploidy is not numeric")}}
      if(length(lambdaPWM)>0){
          if(class(lambdaPWM)=="numeric"){
      NewLambda<-lambdaPWM
      .Object@lambdaPWM<-NewLambda
      validObject(.Object)
      } else {
      stop("lambdaPWM is not numeric
      or a vector of numeric values")
      }}
      if(length(boundMolecules)>0){
          if(class(boundMolecules)=="numeric"){
      NewboundMolecules<-boundMolecules
      .Object@boundMolecules <- NewboundMolecules
      validObject(.Object)
      } else {
      stop("boundMolecules is not numeric")}}

      if(length(backgroundSignal)>0){
          if(class(backgroundSignal)=="numeric"){
      NewbackgroundSignal<-backgroundSignal
      .Object@backgroundSignal <- NewbackgroundSignal
      validObject(.Object)
      } else {
      stop("backgroundSignal is not numeric")}}

      if(length(maxSignal)>0){
          if(class(maxSignal)=="numeric"){
      NewmaxSignal<-maxSignal
      .Object@maxSignal <- NewmaxSignal
      validObject(.Object)
      } else {
      stop("maxSignal is not numeric")}}

      if(length(lociWidth)>0){
          if(class(lociWidth)=="numeric"){
      NewmlociWidth<-lociWidth
      .Object@lociWidth <- NewmlociWidth
      validObject(.Object)
      } else {
      stop("maxSignal is not numeric")}}

      if(length(chipMean)>0){
          if(class(chipMean)=="numeric"){
      NewchipMean<-chipMean
      .Object@chipMean <- NewchipMean
      validObject(.Object)
      } else {
      stop("chipMean is not numeric")}}

      if(length(chipSd)>0){
          if(class(chipSd)=="numeric"){
      NewchipSd<-chipSd
      .Object@chipSd <- NewchipSd
      validObject(.Object)
      } else {
      stop("chipSd is not numeric")}}

      if(length(chipSmooth)>0){
          if(class(chipSmooth)=="numeric"){
      NewchipSmooth<-chipSmooth
      .Object@chipSmooth <- NewchipSmooth
      validObject(.Object)
      }else{
      stop("chipSmooth is not numeric")
      }}
      if(length(stepSize)>0){
          if(class(stepSize)=="numeric"){
      NewstepSize<-stepSize
      .Object@stepSize <- NewstepSize
      validObject(.Object)
      } else {
      stop("stepSize is not numeric")}}

      if(length(removeBackground)>0){
          if(class(removeBackground)=="numeric"){
      NewremoveBackground<-removeBackground
      .Object@removeBackground <- NewremoveBackground
      validObject(.Object)
      } else {
      stop("removeBackground is not numeric")}}
      if(length(PWMpseudocount)>0){
          if(class(PWMpseudocount)=="numeric"){
      Newpseudocount<-PWMpseudocount
      .Object@PWMpseudocount<-Newpseudocount
      validObject(.Object)
      } else {
      stop("PWMpseudocount is not numeric")
      }}

      if(length(naturalLog)>0 & is.logical(naturalLog)){
      NewnaturalLog<-naturalLog
      .Object@naturalLog<-NewnaturalLog
      validObject(.Object)
      }
      if(length(noOfSites)>0){
          if(class(noOfSites)=="numeric"|class(noOfSites)=="character"){
      NewnoOfSites<-noOfSites
      .Object@noOfSites<-NewnoOfSites
      validObject(.Object)
      } else {
      stop("noOfSites is not numeric or 'all' - for all sites")
      }}
      if(length(PWMThreshold)>0){
          if(class(PWMThreshold)=="numeric"){
      NewPWMThreshold<-PWMThreshold
      .Object@PWMThreshold<-NewPWMThreshold
      validObject(.Object)
      } else {
      stop("PWMThreshold is not numeric")
      }}

      if(length(strandRule)>0){
          if(class(strandRule)=="character"){
      NewStrandRule<-strandRule
      .Object@strandRule<-NewStrandRule
      validObject(.Object)
      } else {
      stop("strandRule is not a character string of one
      of the following: max, mean or sum")
      }}

      if(length(whichstrand)>0){
          if(class(whichstrand)=="character"){
      NewStrand<-whichstrand
      .Object@whichstrand<-NewStrand
      validObject(.Object)
      } else {
      stop("whichstrand is not a character string of one of the
      following: +-, -+ , - or +")
      }}
      if(length(noiseFilter)>0){
          if(class(noiseFilter)=="character"){

      .Object@noiseFilter<-noiseFilter
      validObject(.Object)
      } else {
      stop("noiseFilter is not a character string of one of the
      following: zero, mean, median or sigmoid")
      }}

  return(.Object)


    }
)
