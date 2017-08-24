##############################################
#############      S4 Methods    #############
##############################################


## initialize Genomic Profile Paramter and run validity checks for GPP
setMethod(f="initialize",
    signature = "genomicProfileParameters",
    definition=function(.Object,
        PFM=.Object@PFM,
        PWM=.Object@PWM,
        ScalingFactorPWM=.Object@ScalingFactorPWM,
        PWMpseudocount=.Object@PWMpseudocount,
        BPFrequency=.Object@BPFrequency,
        naturalLog=.Object@naturalLog,
        noOfSites=.Object@noOfSites,
        PWMThreshold=.Object@PWMThreshold,
        minPWMScore=.Object@minPWMScore,
        maxPWMScore=.Object@maxPWMScore,
        strandRule=.Object@strandRule,
        whichstrand=.Object@whichstrand){
            if(class(BPFrequency) != "vector" &
            (class(BPFrequency) == "BSgenome" |
            class(BPFrequency) == "DNAStringSet")){
            BPFrequency<-.computeBPFrequency(BPFrequency)
            }
            if(!is.null(BPFrequency) & class(BPFrequency)=="vector") {
            NewBPFrequency<-BPFrequency
            .Object@BPFrequency<-NewBPFrequency
            validObject(.Object)
            }
            if(length(ScalingFactorPWM)>0){
            NewLambda<-ScalingFactorPWM
            .Object@ScalingFactorPWM<-NewLambda
            validObject(.Object)
            }
            if(length(PWMpseudocount)>0){
            Newpseudocount<-PWMpseudocount
            .Object@PWMpseudocount<-Newpseudocount
            validObject(.Object)
            }
            if(length(BPFrequency)>0){
            NewBPFrequency<-BPFrequency
            .Object@BPFrequency<-NewBPFrequency
            validObject(.Object)
            }
            if(length(naturalLog)>0){
            NewnaturalLog<-naturalLog
            .Object@naturalLog<-NewnaturalLog
            validObject(.Object)
            }
            if(length(noOfSites)>0){
            NewnoOfSites<-noOfSites
            .Object@noOfSites<-NewnoOfSites
            validObject(.Object)
            }
            if(length(PWMThreshold)>0){
            NewPWMThreshold<-PWMThreshold
            .Object@PWMThreshold<-NewPWMThreshold
            validObject(.Object)
            }
            if(length(maxPWMScore)>0){
            NewmaxPWMScore<-maxPWMScore
            .Object@maxPWMScore<-NewmaxPWMScore
            validObject(.Object)
            }
            if(length(minPWMScore)>0){
            NewminPWMScore<-minPWMScore
            .Object@minPWMScore<-NewminPWMScore
            validObject(.Object)
            }
            if(class(PFM)!= "matrix" & !is.null(PFM)){
            PFMmatrix<-.parsePFM(PFM)
            .Object@PFM<-PFMmatrix
            validObject(.Object)
            }
            if(any(c("max","mean","sum")==strandRule)){
            NewStrandRule<-strandRule
            .Object@strandRule<-NewStrandRule
            validObject(.Object)
            }
            if(any(c("+","-","+-","-+")==whichstrand)){
            NewStrand<-whichstrand
            .Object@whichstrand<-NewStrand
            validObject(.Object)
            }
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
    return(.Object)
    }
)

setMethod(f="initialize",
    signature = "occupancyProfileParameters",
    definition = function(.Object,
        ploidy = .Object@ploidy,
        boundMolecules = .Object@boundMolecules,
        backgroundSignal = .Object@backgroundSignal,
        maxSignal = .Object@maxSignal,
        chipMean = .Object@chipMean,
        chipSd = .Object@chipSd,
        chipSmooth = .Object@chipSmooth,
        stepSize = .Object@stepSize,
        removeBackground = .Object@removeBackground,
        thetaThreshold = .Object@thetaThreshold){

        if(length(ploidy)>0){
        NewPloidy<-ploidy
        .Object@ploidy <- NewPloidy
        validObject(.Object)
        }
        if(length(boundMolecules)>0){
        NewboundMolecules<-boundMolecules
        .Object@boundMolecules <- NewboundMolecules
        validObject(.Object)
        }
        if(length(backgroundSignal)>0){
        NewbackgroundSignal<-backgroundSignal
        .Object@backgroundSignal <- NewbackgroundSignal
        validObject(.Object)
        }
        if(length(maxSignal)>0){
        NewmaxSignal<-maxSignal
        .Object@maxSignal <- NewmaxSignal
        validObject(.Object)
        }
        if(length(chipMean)>0){
        NewchipMean<-chipMean
        .Object@chipMean <- NewchipMean
        validObject(.Object)
        }
        if(length(chipSd)>0){
        NewchipSd<-chipSd
        .Object@chipSd <- NewchipSd
        validObject(.Object)
        }
        if(length(chipSmooth)>0){
        NewchipSmooth<-chipSmooth
        .Object@chipSmooth <- NewchipSmooth
        validObject(.Object)
        }
        if(length(stepSize)>0){
        NewstepSize<-stepSize
        .Object@stepSize <- NewstepSize
        validObject(.Object)
        }
        if(length(removeBackground)>0){
        NewremoveBackground<-removeBackground
        .Object@removeBackground <- NewremoveBackground
        validObject(.Object)
        }
        if(length(thetaThreshold)>0){
        NewthetaThreshold<-thetaThreshold
        .Object@thetaThreshold <- NewthetaThreshold
        validObject(.Object)
        }
    return(.Object)
    }
)

### initialize Genomic Profile by generating a PWM from a PFM
setMethod(".generatePWM","genomicProfileParameters",
    function(object){
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
