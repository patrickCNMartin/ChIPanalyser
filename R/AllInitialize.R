##############################################
#############      S4 Methods    #############
##############################################


## initialize Genomic Profile Paramter and run validity checks for GPP
setMethod(f="initialize",
    signature = "genomicProfileParameters",
    definition=function(.Object,
        PFM=.Object@PFM,
        PWM=.Object@PWM,
        PFMFormat = .Object@PFMFormat,
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
                if(class(ScalingFactorPWM)=="numeric"){
            NewLambda<-ScalingFactorPWM
            .Object@ScalingFactorPWM<-NewLambda
            validObject(.Object)
            } else {
            stop("ScalingFactorPWM is not numeric
            or a vector of numeric values")
            }}
            if(length(PFMFormat)>0){
                if(class(PFMFormat)=="character"){
            NewPFMFormat<-PFMFormat
            .Object@PFMFormat<-NewPFMFormat
            validObject(.Object)
            } else {
            stop("PFMFormat is not a character string of one of the following:
            raw, transfac, JASPAR, sequences or matrix")
            }}
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
                if(class(noOfSites)=="numeric"){
            NewnoOfSites<-noOfSites
            .Object@noOfSites<-NewnoOfSites
            validObject(.Object)
            } else {
            stop("noOfSites is not numeric")
            }}
            if(length(PWMThreshold)>0){
                if(class(PWMThreshold)=="numeric"){
            NewPWMThreshold<-PWMThreshold
            .Object@PWMThreshold<-NewPWMThreshold
            validObject(.Object)
            } else {
            stop("PWMThreshold is not numeric")
            }}
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
            if(class(ploidy)=="numeric"){
        NewPloidy<-ploidy
        .Object@ploidy <- NewPloidy
        validObject(.Object)
        } else {
        stop("ploidy is not numeric")}}

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

        if(length(thetaThreshold)>0){
            if(class(thetaThreshold)=="numeric"){
        NewthetaThreshold<-thetaThreshold
        .Object@thetaThreshold <- NewthetaThreshold
        validObject(.Object)
        } else {
        stop("thetaThreshold is not numeric")}}

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
